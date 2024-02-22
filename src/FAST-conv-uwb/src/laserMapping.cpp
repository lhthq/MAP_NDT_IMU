// This is an advanced implementation of the algorithm described in the
// following paper:
//   J. Zhang and S. Singh. LOAM: Lidar Odometry and Mapping in Real-time.avr
//     Robotics: Science and Systems Conference (RSS). Berkeley, CA, July 2014.

// Modifier: Livox               dev@livoxtech.com

// Copyright 2013, Ji Zhang, Carnegie Mellon University
// Further contributions copyright (c) 2016, Southwest Research Institute
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
#include "udp_comm.h"
#include <Eigen/Core>
#include <Python.h>
#include <csignal>
#include <fstream>
#include <math.h>
#include <mutex>
#include <omp.h>
#include <ros/ros.h>
#include <so3_math.h>
#include <thread>
#include <unistd.h>
// #include <common_lib.h>
#include "IMU_Processing.h"
#include "laserMapping.h"
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>

#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <sensor_msgs/PointCloud2.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#include <pcl/registration/ndt.h>

#include "preprocess.h"

#include <geometry_msgs/Vector3.h>
#include <livox_ros_driver/CustomMsg.h>
#include <opencv2/opencv.hpp>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <tf2/transform_datatypes.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2_eigen/tf2_eigen.h>
#include <tf2_ros/transform_listener.h>

#include "color.h"

#include <sensor_msgs/NavSatFix.h>

#include "ndt_3d.h"

std::vector<Eigen::Vector3d> lla_vec;

pcl::NormalDistributionsTransform<pcl::PointXYZ, pcl::PointXYZ> pclndt_;

PointCloudXYZI::Ptr merge_points(new PointCloudXYZI());
pcl::PointCloud<pcl::PointXYZ>::Ptr map_points_ptr(new pcl::PointCloud<pcl::PointXYZ>);
double record_time = 0;
double last_record_time = 0;

bool have_gnd_normal = false;

bool lio_work = false;

double map_adjust_rate = 0.1;
StatesGroup state_lio;

#ifdef USE_ikdtree
#ifdef USE_ikdforest
#include <ikd-Forest/ikd_Forest.h>
#else
#include <ikd-Tree/ikd_Tree.h>
#endif
#else
#include <pcl/kdtree/kdtree_flann.h>
#endif

#define INIT_TIME (0.0)
#define MAXN (360000)
#define PUBFRAME_PERIOD (20)
#define WINDOW_SIZE (25)

float DET_RANGE = 300.0f;
#ifdef USE_ikdforest
const int laserCloudWidth = 200;
const int laserCloudHeight = 200;
const int laserCloudDepth = 200;
const int laserCloudNum = laserCloudWidth * laserCloudHeight * laserCloudDepth;

#else
const float MOV_THRESHOLD = 1.5f;
#endif

mutex mtx_buffer;
condition_variable sig_buffer;

Eigen::Vector3d zero_utm;
Eigen::Vector3d zero_utm_imu;
// V3D init_enu;

Eigen::Vector3d zero_utm_pub;
bool gps_init = true;
uint8_t zone_for_pub;
char band_for_pub;
float get_gps_conv_X;
float get_gps_conv_Y;
float get_gps_conv_Z;
double last_timestamp_gnss;
bool first_init = true;
ros::Publisher pubGNSSIMUGT;
std::deque<sensor_msgs::NavSatFix> gps_queue;


std::deque<nav_msgs::Odometry> gtQueue;


bool init_map_wx = false;

std::mutex ndt_map_mtx_;
Eigen::Matrix4f pre_trans, delta_trans, initial_pose_matrix;

int turn = 0;
int lidarswitch = 0;
int gnssswitch = 0;

Eigen::Affine3f affine_lidar_to_imu;
Eigen::Affine3f affine_imu_to_body;
Eigen::Affine3f affine_lidar_to_body;

Eigen::Affine3f affine_imu_to_map;


double t_lidar = 0.0;
double t_uwb = 0.0;

double t_lastnhc = 0.0;
double t_last_frame = 0.0;

double yaw_dual_rtk = 0.0;
int avr_sign = 0;

string root_dir = ROOT_DIR;
string map_file_path, lid_topic, imu_topic, hilti_seq_name, img_topic, gnss_topic,
    config_file;

M3D Eye3d(M3D::Identity());
M3F Eye3f(M3F::Identity());
V3D Zero3d(0, 0, 0);
V3F Zero3f(0, 0, 0);
V3D extT(Zero3d);
V3D gilever(Zero3d);
V3D zupt_std(Zero3d);
Vector2d nhc_std(0.0, 0.0);
M3D extR(Eye3d);
M3D _gravity_correct_rotM(M3D::Identity());


vector<Pose6D> imu_vector;
vector<imu_range> imu_ra;

int effctive_u = 0;

int iterCount = 0, feats_down_size = 0, NUM_MAX_ITERATIONS = 0,
    laserCloudValidNum = 0, effct_feat_num = 0, time_log_counter = 0,
    publish_count = 0;
int MIN_IMG_COUNT = 0;

double res_mean_last = 0.05, _last_lidar_processed_time = -1.0;
double gyr_cov = 0, acc_cov = 0;
double last_timestamp_lidar = -1.0, last_timestamp_imu = -1.0,
       last_timestamp_img = -1.0;
double filter_size_corner_min = 0, filter_size_surf_min = 0,
       filter_size_map_min = 0, fov_deg = 0;
double cube_len = 0, HALF_FOV_COS = 0, FOV_DEG = 0, total_distance = 0,
       lidar_end_time = 0, first_lidar_time = 0.0;
double first_img_time = -1.0;
double kdtree_incremental_time = 0, kdtree_search_time = 0,
       kdtree_delete_time = 0.0;
int kdtree_search_counter = 0, kdtree_size_st = 0, kdtree_size_end = 0,
    add_point_size = 0, kdtree_delete_counter = 0;
;
double copy_time = 0, readd_time = 0, fov_check_time = 0, readd_box_time = 0,
       delete_box_time = 0;
double T1[MAXN], T2[MAXN], s_plot[MAXN], s_plot2[MAXN], s_plot3[MAXN],
    s_plot4[MAXN], s_plot5[MAXN], s_plot6[MAXN], s_plot7[MAXN];
double keyf_rotd = 0.05, keyf_posd = 0.15;
double match_time = 0, solve_time = 0, solve_const_H_time = 0;

/*** For voxel map ***/
double max_voxel_size, min_eigen_value = 0.003, match_s = 0.90, sigma_num = 2.0,
                       match_eigen_value = 0.0025;
double beam_err = 0.03, dept_err = 0.05;
int pub_map = 0, voxel_layer = 1;
int last_match_num = 0;
bool init_map = false, use_new_map = true, is_pub_plane_map = false,
     pcd_save_en = false, img_save_en = false, USE_NED = true,
     effect_point_pub = false, hilti_en = false;
int min_points_size = 30, pcd_save_type = 0, pcd_save_interval = -1,
    img_save_interval = 1, pcd_index = 0;
std::time_t startTime, endTime;
std::unordered_map<VOXEL_LOC, OctoTree *> feat_map;
V3D layer_size(20, 10, 10);
std::vector<M3D> crossmat_list;


bool lidar_pushed, imu_en, flg_reset, flg_exit = false;
int dense_map_en = 1;
int img_en = 1, imu_int_frame = 3;
int lidar_en = 1;
int uwb_en = 1;
int gnss_en = 1;
int debug = 0;
bool is_first_frame = false;
int grid_size, patch_size;
double outlier_threshold;
int udb_comm_en = 0;
std::string udp_ip = "BROADCAST";
std::string udp_port = "8000";

vector<BoxPointType> cub_needrm;
vector<BoxPointType> cub_needad;

// deque<sensor_msgs::PointCloud2::ConstPtr> lidar_buffer;
deque<PointCloudXYZI::Ptr> lidar_buffer;
deque<double> time_buffer;
deque<sensor_msgs::Imu::ConstPtr> imu_buffer;
deque<cv::Mat> img_buffer;
deque<double> img_time_buffer;
deque<double> uwb_time_buffer;


vector<bool> point_selected_surf;
vector<vector<int>> pointSearchInd_surf;
vector<PointVector> Nearest_Points;
vector<double> res_last;
vector<double> extrinT(3, 0.0);
vector<double> ginslever(3, 0.0);


vector<double> extrinR(9, 0.0);
vector<double> cameraextrinT(3, 0.0);
vector<double> cameraextrinR(9, 0.0);
double total_residual;
double LASER_POINT_COV, IMG_POINT_COV, cam_fx, cam_fy, cam_cx, cam_cy;
bool flg_EKF_inited, flg_EKF_converged, EKF_stop_flg = 0;
// surf feature in map
PointCloudXYZI::Ptr featsFromMap(new PointCloudXYZI());
PointCloudXYZI::Ptr cube_points_add(new PointCloudXYZI());
PointCloudXYZI::Ptr map_cur_frame_point(new PointCloudXYZI());
PointCloudXYZI::Ptr sub_map_cur_frame_point(new PointCloudXYZI());

PointCloudXYZI::Ptr feats_undistort(new PointCloudXYZI());
PointCloudXYZI::Ptr feats_down_body(new PointCloudXYZI());
PointCloudXYZI::Ptr feats_down_world(new PointCloudXYZI());
PointCloudXYZI::Ptr normvec(new PointCloudXYZI(100000, 1));
PointCloudXYZI::Ptr laserCloudOri(new PointCloudXYZI(100000, 1));
PointCloudXYZI::Ptr corr_normvect(new PointCloudXYZI(100000, 1));

ofstream fout_pre, fout_out, fout_dbg, fout_pcd_pos, fout_img_pos;
ofstream evo_file;
FILE *fp_rtk;
string log_dirrtk = root_dir + "/Log/rtk.txt";


pcl::VoxelGrid<PointType> downSizeFilterSurf;
pcl::VoxelGrid<PointType> downSize;

#ifdef USE_ikdtree
#ifdef USE_ikdforest
KD_FOREST ikdforest;
#else
KD_TREE ikdtree;
#endif
#else
pcl::KdTreeFLANN<PointType>::Ptr
    kdtreeSurfFromMap(new pcl::KdTreeFLANN<PointType>());
#endif

V3F XAxisPoint_body(LIDAR_SP_LEN, 0.0, 0.0);
V3F XAxisPoint_world(LIDAR_SP_LEN, 0.0, 0.0);
V3D euler_cur;
V3D position_last(Zero3d);
V3D position_last_ndt(Zero3d);
Eigen::Matrix3d Rcl;
Eigen::Vector3d Pcl;

LidarMeasureGroup LidarMeasures;

#ifdef USE_IKFOM
esekfom::esekf<state_ikfom, 12, input_ikfom> kf;
state_ikfom state_point;
vect3 pos_lid;
#else
StatesGroup state;
StatesGroup last_u_state;
#endif

nav_msgs::Path path;
nav_msgs::Path rtkpath;

nav_msgs::Odometry odomAftMapped;
geometry_msgs::Quaternion geoQuat;
geometry_msgs::PoseStamped msg_body_pose;

shared_ptr<Preprocess> p_pre(new Preprocess());

void SigHandle(int sig)
{
  flg_exit = true;
  ROS_WARN("catch sig %d", sig);
  sig_buffer.notify_all();
}

void pointBodyToWorld(const PointType &pi, PointType &po)
{
  V3D p_body(pi.x, pi.y, pi.z);

  V3D p_global(state.rot_end * (extR * p_body + extT) + state.pos_end);


  po.x = p_global(0);
  po.y = p_global(1);
  po.z = p_global(2);
  po.intensity = pi.intensity;
}

template <typename T>
void pointBodyToWorld(const Matrix<T, 3, 1> &pi, Matrix<T, 3, 1> &po)
{
  V3D p_body(pi[0], pi[1], pi[2]);

  V3D p_global(state.rot_end * (extR * p_body + extT) + state.pos_end);

  po[0] = p_global(0);
  po[1] = p_global(1);
  po[2] = p_global(2);
}

template <typename T>
Matrix<T, 3, 1> pointBodyToWorld(const Matrix<T, 3, 1> &pi)
{
  V3D p(pi[0], pi[1], pi[2]);

  p = (state.rot_end * (extR * p + extT) + state.pos_end);

  Matrix<T, 3, 1> po(p[0], p[1], p[2]);

  return po;
}

void frameBodyToWorld(const PointCloudXYZI::Ptr &pi, PointCloudXYZI::Ptr &po)
{
  int pi_size = pi->points.size();
  po->resize(pi_size);
  for (int i = 0; i < pi_size; i++)
  {
    /* transform to world frame */
    pointBodyToWorld(pi->points[i], po->points[i]);
  }
}

void RGBpointBodyToWorld(PointType const *const pi, PointType *const po)
{
  V3D p_body(pi->x, pi->y, pi->z);

  V3D p_global(state.rot_end * (extR * p_body + extT) + state.pos_end);


  p_global = _gravity_correct_rotM * p_global;

  po->x = p_global(0);
  po->y = p_global(1);
  po->z = p_global(2);
  po->intensity = pi->intensity;

  float intensity = pi->intensity;
  intensity = intensity - floor(intensity);

  int reflection_map = intensity * 10000;
}

void RGBpointBodyLidarToIMU(PointType const *const pi, PointType *const po)
{
  V3D p_body_lidar(pi->x, pi->y, pi->z);

  V3D p_body_imu(extR * p_body_lidar + extT);


  po->x = p_body_imu(0);
  po->y = p_body_imu(1);
  po->z = p_body_imu(2);
  po->intensity = pi->intensity;
}

void standard_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  if (!lidar_en)
    return;
  mtx_buffer.lock();
  //  cout<<"got feature"<<endl;
  if (msg->header.stamp.toSec() < last_timestamp_lidar)
  {
    ROS_ERROR("lidar loop back, clear buffer");
    lidar_buffer.clear();
  }
  // ROS_INFO("get point cloud at time: %.6f", msg->header.stamp.toSec());
  PointCloudXYZI::Ptr ptr(new PointCloudXYZI());
  p_pre->process(msg, ptr);
  lidar_buffer.push_back(ptr);
  time_buffer.push_back(msg->header.stamp.toSec());
  last_timestamp_lidar = msg->header.stamp.toSec();

  mtx_buffer.unlock();
  sig_buffer.notify_all();
}

void livox_pcl_cbk(const livox_ros_driver::CustomMsg::ConstPtr &msg)
{
  if (!lidar_en)
    return;
  mtx_buffer.lock();
  ROS_INFO("get LiDAR, its header time: %.6f", msg->header.stamp.toSec());
  if (msg->header.stamp.toSec() < last_timestamp_lidar)
  {
    ROS_ERROR("lidar loop back, clear buffer");
    lidar_buffer.clear();
  }
  // ROS_INFO("get point cloud at time: %.6f", msg->header.stamp.toSec());
  PointCloudXYZI::Ptr ptr(new PointCloudXYZI());
  p_pre->process(msg, ptr);
  lidar_buffer.push_back(ptr);
  time_buffer.push_back(msg->header.stamp.toSec());
  last_timestamp_lidar = msg->header.stamp.toSec();

  mtx_buffer.unlock();
  sig_buffer.notify_all();
}


void imu_cbk(const sensor_msgs::Imu::ConstPtr &msg_in)
{
  if (!imu_en)
    return;

  if (last_timestamp_lidar < 0.0)
    return;
  publish_count++;
  // ROS_INFO("get imu at time: %.6f", msg_in->header.stamp.toSec());
  sensor_msgs::Imu::Ptr msg(new sensor_msgs::Imu(*msg_in));

  double timestamp = msg->header.stamp.toSec();
  mtx_buffer.lock();

  if (last_timestamp_imu > 0.0 && timestamp < last_timestamp_imu)
  {
    mtx_buffer.unlock();
    sig_buffer.notify_all();
    ROS_ERROR("imu loop back \n");
    return;
  }
  // old 0.2
  if (last_timestamp_imu > 0.0 && timestamp > last_timestamp_imu + 0.4)
  {
    mtx_buffer.unlock();
    sig_buffer.notify_all();
    ROS_WARN("imu time stamp Jumps %0.4lf seconds \n",
             timestamp - last_timestamp_imu);
    return;
  }

  last_timestamp_imu = timestamp;

  imu_buffer.push_back(msg);
  // cout<<"got imu: "<<timestamp<<" imu size "<<imu_buffer.size()<<endl;
  mtx_buffer.unlock();
  sig_buffer.notify_all();
}

void gtHandler(const nav_msgs::Odometry::ConstPtr &gtMsg)
{
  if (isnan(gtMsg->pose.pose.position.x) || isnan(gtMsg->pose.pose.position.y) || isnan(gtMsg->pose.pose.position.z))
  {
    return;
  }
  // cout<<"ground truth recieve"<<endl;
  Eigen::Affine3f trans_body_to_map, trans_imu_to_map;
  odometryMsgToAffine3f(*gtMsg, trans_body_to_map);
  trans_imu_to_map = trans_body_to_map * affine_imu_to_body;
  float roll, pitch, yaw, x, y, z;
  pcl::getTranslationAndEulerAngles(trans_imu_to_map, x, y, z, roll, pitch, yaw);
  tf::Quaternion q = tf::createQuaternionFromRPY(roll, pitch, yaw);
  nav_msgs::Odometry gtMsgNew;
  gtMsgNew.header = gtMsg->header;
  gtMsgNew.pose.covariance = gtMsg->pose.covariance;
  gtMsgNew.pose.pose.position.x = x;
  gtMsgNew.pose.pose.position.y = y;
  gtMsgNew.pose.pose.position.z = z;
  gtMsgNew.pose.pose.orientation.x = q.x();
  gtMsgNew.pose.pose.orientation.y = q.y();
  gtMsgNew.pose.pose.orientation.z = q.z();
  gtMsgNew.pose.pose.orientation.w = q.w();
  gtQueue.push_back(gtMsgNew);

}



double img_lid_timediff = 0;


//  without LiDAR scan ->leverage UWB/INS integration  liuhong20230414
bool sync_packages(LidarMeasureGroup &meas)
{
  if ((lidar_buffer.empty() && lidar_en) || (img_buffer.empty() && img_en) || (gps_queue.empty() && gnss_en))
  {
    // cout<<"lidar_en:"<<lidar_en<<"   "<<"gnss_en:"<<gnss_en<<endl;
    return false;
  }
  if (meas.is_lidar_end) // If meas.is_lidar_end==true, means it just after scan    //first scan is false
                         // end, clear all buffer in meas.
  {
    meas.measures.clear();
    meas.is_lidar_end = false;
  }

  if (!lidar_pushed) // first scan is false
  {                  // If not in lidar scan, need to generate new meas
    if (lidar_buffer.empty())
    {
      // ROS_ERROR("out sync");
      return false;
    }
    meas.lidar = lidar_buffer.front(); // push the firsrt lidar topic
    if (meas.lidar->points.size() <= 1)
    {
      mtx_buffer.lock();

      if (gps_queue.size() > 0) // temp method, ignore uwb topic when no lidar points, keep sync
      {
        lidar_buffer.pop_front();
        gps_queue.pop_front();
      }
      mtx_buffer.unlock();
      sig_buffer.notify_all();
      // ROS_ERROR("out sync");
      return false;
    }
    // sort(meas.lidar->points.begin(), meas.lidar->points.end(), time_list); //
    // sort by sample timestamp
    meas.lidar_beg_time = time_buffer.front();                                                 // generate lidar_beg_time
    lidar_end_time = meas.lidar_beg_time + meas.lidar->points.back().curvature / double(1000); // calc lidar scan end time
    lidar_pushed = true;                                                                       // flag

    // if (!imu_en && !img_en) // without imu and uwb
    if (!imu_en && !gnss_en) // no imu and uwb topic, means   LO  system
    {
      cout << "got points, imu and image disabled!" << endl;
      return true;
    }
  }

  // if (img_buffer.empty())
  // if (uwb_buffer.empty()) // no uwb topic, means   LIO  system  取当前LiDAR最近的range信息，会浪费中间UWB的量测数据（50HZ）
  {
    if (imu_en && last_timestamp_imu < lidar_end_time) // last_timestamp_imu  最新的IMU
    {                                                  // imu message needs to be larger than
      // lidar_end_time, keep complete propagate.
      // ROS_ERROR("out sync");
      return false;
    }

    struct MeasureGroup m; // standard method to keep imu message.

    if (!imu_buffer.empty())
    {
      double imu_time = imu_buffer.front()->header.stamp.toSec();
      m.imu.clear();
      mtx_buffer.lock();

      while ((!imu_buffer.empty() && (imu_time < lidar_end_time)))
      {
        imu_time = imu_buffer.front()->header.stamp.toSec();
        if (imu_time > lidar_end_time)
          break;
        m.imu.push_back(imu_buffer.front());
        imu_buffer.pop_front();
      }
    }

    if (!gps_queue.empty())
    {
      double time_g = gps_queue.front().header.stamp.toSec();
      m.range_in.clear();
      mtx_buffer.lock();
      while ((!gps_queue.empty() && (time_g < lidar_end_time)))
      {
        time_g = gps_queue.front().header.stamp.toSec();
        if (time_g > lidar_end_time)
          break;
        // if (uwb_buffer.front().id == 0) // node 0
        //{
        m.range_in.push_back(gps_queue.front());
        //}

        gps_queue.pop_front();
      }
    }

    lidar_buffer.pop_front();
    time_buffer.pop_front();
    mtx_buffer.unlock();

    sig_buffer.notify_all();
    lidar_pushed = false;     // sync one whole lidar scan.
    meas.is_lidar_end = true; // process lidar topic, so timestamp should be lidar scan end.
    meas.measures.push_back(m);
    // ROS_INFO("ONlY HAS LiDAR and IMU, NO IMAGE!");
    cout << "imu_buffer.size(): " << m.imu.size() << " gnss_buffer.size(): " << m.range_in.size() << endl;
    _last_lidar_processed_time = meas.lidar_beg_time + meas.lidar->points.back().curvature / double(1000);
    return true;
  }
}


#define NUM_POINTS 2000

// PointCloudXYZRGB::Ptr pcl_wait_pub_RGB(new PointCloudXYZRGB(500000, 1));
PointCloudXYZI::Ptr pcl_wait_pub(new PointCloudXYZI(500000, 1));
PointCloudXYZI::Ptr pcl_wait_save(new PointCloudXYZI());

StatesGroup last_state;


template <typename T>
void set_posestamp(T &out)
{
#ifdef USE_IKFOM
  // state_ikfom stamp_state = kf.get_x();
  out.position.x = state_point.pos(0);
  out.position.y = state_point.pos(1);
  out.position.z = state_point.pos(2);
#else
  out.position.x = state.pos_end(0);
  out.position.y = state.pos_end(1);
  out.position.z = state.pos_end(2);
#endif
  out.orientation.x = geoQuat.x;
  out.orientation.y = geoQuat.y;
  out.orientation.z = geoQuat.z;
  out.orientation.w = geoQuat.w;
}

void publish_odometry(const ros::Publisher &pubOdomAftMapped)
{
  odomAftMapped.header.frame_id = "camera_init";
  odomAftMapped.child_frame_id = "aft_mapped";
  //odomAftMapped.header.stamp = ros::Time()fromSec(last_timestamp_lidar);
  odomAftMapped.header.stamp =  ros::Time().fromSec(LidarMeasures.lidar_beg_time+0.1);

   //   ros::Time::now(); //.ros::Time()fromSec(last_timestamp_lidar);
  set_posestamp(odomAftMapped.pose.pose);

  static tf::TransformBroadcaster br;
  tf::Transform transform;
  tf::Quaternion q;
  transform.setOrigin(tf::Vector3(state.pos_end(0), state.pos_end(1),
                                  state.pos_end(2)));
  q.setW(geoQuat.w);
  q.setX(geoQuat.x);
  q.setY(geoQuat.y);
  q.setZ(geoQuat.z);
  transform.setRotation(q);
  br.sendTransform(tf::StampedTransform(transform,
                                        odomAftMapped.header.stamp, "camera_init", "aft_mapped"));
  pubOdomAftMapped.publish(odomAftMapped);
}


void publish_path(const ros::Publisher pubPath)
{
  set_posestamp(msg_body_pose.pose);
  msg_body_pose.header.stamp = ros::Time::now();
  msg_body_pose.header.frame_id = "camera_init";
  static int jjj = 0;
  jjj++;
  if (jjj % 1 == 0)
  {
    path.poses.push_back(msg_body_pose);
    pubPath.publish(path);
  }
}

void readParameters(ros::NodeHandle &nh)
{
  nh.param<int>("dense_map_enable", dense_map_en, 1);
  nh.param<int>("img_enable", img_en, 1);
  nh.param<int>("uwb_enable", uwb_en, 1);

  nh.param<int>("gnss_enable", gnss_en, 1);

  nh.param<int>("lidar_enable", lidar_en, 1);
  nh.param<int>("udb_comm_enable", udb_comm_en, 0);
  nh.param<int>("debug", debug, 0);
  nh.param<int>("max_iteration", NUM_MAX_ITERATIONS, 4);
  nh.param<int>("min_img_count", MIN_IMG_COUNT, 1000);

  nh.param<double>("cam_fx", cam_fx, 453.483063);
  nh.param<double>("cam_fy", cam_fy, 453.254913);
  nh.param<double>("cam_cx", cam_cx, 318.908851);
  nh.param<double>("cam_cy", cam_cy, 234.238189);

  nh.param<double>("laser_point_cov", LASER_POINT_COV, 0.001);
  nh.param<double>("img_point_cov", IMG_POINT_COV, 10);
  nh.param<string>("map_file_path", map_file_path, "");
  nh.param<string>("udp_ip", udp_ip, "BROADCAST");
  nh.param<string>("udp_port", udp_port, "8000");
  nh.param<string>("common/lid_topic", lid_topic, "/livox/lidar");
  nh.param<string>("common/imu_topic", imu_topic, "/livox/imu");
  nh.param<string>("hilti/seq", hilti_seq_name, "01");
  nh.param<bool>("hilti/en", hilti_en, false);
  nh.param<string>("camera/img_topic", img_topic, "/usb_cam/image_raw");
  nh.param<double>("filter_size_corner", filter_size_corner_min, 0.5);
  nh.param<double>("filter_size_surf", filter_size_surf_min, 0.5);
  nh.param<double>("filter_size_map", filter_size_map_min, 0.5);
  nh.param<double>("cube_side_length", cube_len, 200);
  nh.param<double>("mapping/fov_degree", fov_deg, 180);
  nh.param<double>("mapping/gyr_cov", gyr_cov, 1.0);
  nh.param<double>("mapping/acc_cov", acc_cov, 1.0);
  nh.param<int>("mapping/imu_int_frame", imu_int_frame, 3);
  nh.param<bool>("mapping/imu_en", imu_en, false);
  nh.param<float>("mapping/det_range", DET_RANGE, 300.f);
  nh.param<bool>("voxel/voxel_map_en", use_new_map, false);
  nh.param<bool>("voxel/pub_plane_en", is_pub_plane_map, false);
  nh.param<double>("voxel/match_eigen_value", match_eigen_value, 0.0025);
  nh.param<int>("voxel/layer", voxel_layer, 1);
  nh.param<double>("voxel/match_s", match_s, 0.90);
  nh.param<double>("voxel/voxel_size", max_voxel_size, 1.0);
  nh.param<double>("voxel/min_eigen_value", min_eigen_value, 0.01);
  nh.param<double>("voxel/sigma_num", sigma_num, 3);
  nh.param<double>("voxel/beam_err", beam_err, 0.02);
  nh.param<double>("voxel/dept_err", dept_err, 0.05);
  nh.param<double>("preprocess/blind", p_pre->blind, 0.01);
  nh.param<double>("image_save/rot_dist", keyf_rotd, 0.01);
  nh.param<double>("image_save/pos_dist", keyf_posd, 0.01);
  nh.param<int>("preprocess/lidar_type", p_pre->lidar_type, AVIA);
  nh.param<int>("preprocess/scan_line", p_pre->N_SCANS, 16);
  nh.param<int>("point_filter_num", p_pre->point_filter_num, 2);
  nh.param<int>("pcd_save/interval", pcd_save_interval, -1);
  nh.param<int>("image_save/interval", img_save_interval, 1);
  nh.param<int>("pcd_save/type", pcd_save_type, 0);
  nh.param<bool>("pcd_save/pcd_save_en", pcd_save_en, false);
  nh.param<bool>("image_save/img_save_en", img_save_en, false);
  nh.param<bool>("feature_extract_enable", p_pre->feature_enabled, false);
  nh.param<vector<double>>("mapping/extrinsic_T", extrinT, vector<double>());
  nh.param<vector<double>>("mapping/lever", ginslever, vector<double>());


  nh.param<vector<double>>("mapping/extrinsic_R", extrinR, vector<double>());
  nh.param<vector<double>>("camera/Pcl", cameraextrinT, vector<double>());
  nh.param<vector<double>>("camera/Rcl", cameraextrinR, vector<double>());
  nh.param<int>("grid_size", grid_size, 40);
  nh.param<int>("patch_size", patch_size, 4);
  nh.param<double>("outlier_threshold", outlier_threshold, 100);
  nh.param<bool>("publish/effect_point_pub", effect_point_pub, false);
  nh.param<string>("common/gnss_topic", gnss_topic, "/gps/fix");

  p_pre->blind_sqr = p_pre->blind * p_pre->blind;

  cout << "ranging cov:" << dept_err << " , angle cov:" << beam_err << "DET_RANGE:" << DET_RANGE
       << std::endl;
}



void loose_observe(StatesGroup &state_propagat, V3D &observe, M3D &observe_rot, double trans_noise, double rot_noise)
{
  V3D p_gnss;
  p_gnss << observe(0), observe(1), observe(2); // 观测

  Vec6d dz;
  M3D rot_diff = (state.rot_end.inverse()) * observe_rot;
  dz.head<3>() = Log(rot_diff);
  dz.tail<3>() = p_gnss - state.pos_end;

  Eigen::MatrixXd R_gnsspos; // 协方差
  Vec6d std;
  std << rot_noise * rot_noise, rot_noise * rot_noise, rot_noise * rot_noise, trans_noise * trans_noise, trans_noise * trans_noise, trans_noise * trans_noise;
  R_gnsspos = std.asDiagonal();

  Eigen::MatrixXd I;
  I.resizeLike(state.cov);
  Eigen::MatrixXd G_gnss;
  G_gnss.resizeLike(state.cov);
  Eigen::MatrixXd K;
  // EKF_stop_flg = false;
  // flg_EKF_converged = false;

  // int iterations = 10;
  // for (int iterCount = 0; iterCount < iterations; iterCount++)
  // {
  // M3D jaco;
  // V3D TEMP_ = state.rot_end * gilever;
  // jaco << SKEW_SYM_MATRX(TEMP_);
  Eigen::MatrixXd H_gnsspos;
  H_gnsspos.resize(6, DIM_STATE); // 2.观测矩阵
  H_gnsspos.setZero();
  H_gnsspos.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
  H_gnsspos.block(3, 3, 3, 3) = Eigen::Matrix3d::Identity();

  auto temp = H_gnsspos * state.cov * H_gnsspos.transpose() + R_gnsspos;
  K = state.cov * H_gnsspos.transpose() * temp.inverse();
  I.setIdentity();
  G_gnss = I - K * H_gnsspos;

  // auto vec = state_propagat - state;

  VD(DIM_STATE)
  sta_temp = K * dz; // + G_gnss * vec
                     // state = state + solution;

  V3D ro(sta_temp(0), sta_temp(1), sta_temp(2));
  M3D deltarom;
  deltarom << SKEW_SYM_MATRX(ro);
  M3D ind;
  ind.setIdentity();

  state = state + (sta_temp);

  // state.rot_end = (ind + deltarom) * state.rot_end;
  // state.pos_end = state.pos_end - sta_temp.block<3, 1>(3, 0);
  // state.vel_end = state.vel_end - sta_temp.block<3, 1>(6, 0);
  // state.bias_g = state.bias_g + sta_temp.block<3, 1>(9, 0);
  // state.bias_a = state.bias_a + sta_temp.block<3, 1>(12, 0);
  // state.gravity = state.gravity + sta_temp.block<3, 1>(15, 0);

  state.cov = G_gnss * state.cov; //* G_gnss.transpose() + K * R_gnsspos * K.transpose()

  //   auto rot_add = (solution).block<3, 1>(0, 0);
  //   auto t_add = (solution).block<3, 1>(3, 0);
  //   printf(YELLOW "[check map rot_add  and t_add] ");
  //   cout << "ROT_add.norm() * 57.3:  " << rot_add.norm() * 57.3 << "   T_add.norm() * 100:  " << t_add.norm() * 100 << endl;
  //   cout << RESET;
  //   if ((rot_add.norm() * 57.3 < 0.01) && (t_add.norm() * 100 < 0.015))
  //   {
  //     flg_EKF_converged = true;
  //   }
  //   if (iterCount == iterations - 1 || flg_EKF_converged == true)
  //   {
  //     EKF_stop_flg = true;
  //     /*** Covariance Update ***/
  //     state.cov = G_gnss * state.cov;
  //     cout << "state  pos:" << state.pos_end << endl;
  //   }
  //   if (EKF_stop_flg)
  //     break;
  // }
}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "laserMapping");
  ros::NodeHandle nh;
  readParameters(nh);
  cout << "debug:" << debug << " MIN_IMG_COUNT: " << MIN_IMG_COUNT << endl;
  pcl_wait_pub->clear();

  ros::Subscriber sub_pcl =
      p_pre->lidar_type == AVIA
          ? nh.subscribe(lid_topic, 200000, livox_pcl_cbk)
          : nh.subscribe(lid_topic, 1, standard_pcl_cbk);


  ros::Subscriber sub_imu = nh.subscribe(imu_topic, 200000, imu_cbk);
  ros::Subscriber sub_gt = nh.subscribe("/ground_truth", 10, gtHandler);

  ros::Publisher pubLaserCloudFullRes =
      nh.advertise<sensor_msgs::PointCloud2>("/cloud_registered", 100);
  ros::Publisher pubOdomAftMapped =
      nh.advertise<nav_msgs::Odometry>("/aft_mapped_to_init", 10);
  ros::Publisher pubPath = nh.advertise<nav_msgs::Path>("/path", 10);


  ros::Publisher plane_pub =
      nh.advertise<visualization_msgs::Marker>("/planner_normal", 1);
  ros::Publisher voxel_pub =
      nh.advertise<visualization_msgs::MarkerArray>("/voxels", 1);

  pubGNSSIMUGT = nh.advertise<nav_msgs::Odometry>("/gnss_imu_gt", 1);

  ros::Publisher sensor_aligned_pose_pub_ = nh.advertise<sensor_msgs::PointCloud2>("points_aligned", 10);


  path.header.stamp = ros::Time::now();
  path.header.frame_id = "camera_init";

  rtkpath.header.stamp = ros::Time::now();
  rtkpath.header.frame_id = "camera_init";

  affine_imu_to_body = pcl::getTransformation(-0.11, -0.18, -0.71, 0.0, 0.0, 0.0);
  affine_lidar_to_body = pcl::getTransformation(0.002, -0.004, -0.957, 0.014084807063594, 0.002897246558311, -1.583065991436417);
  affine_lidar_to_imu = affine_imu_to_body.inverse() * affine_lidar_to_body;


  string saveDirectory = root_dir + "/Log/";
  string sequence = "KML";
  string configDirectory = root_dir + "/config/";

  std::string path = "/home/seu/c7map_for_wx/2012-01-15/map_pcd/GlobalMap.pcd";
  sensor_msgs::PointCloud2 pcd;
  pcl::io::loadPCDFile(path.c_str(), pcd);
  pcl::fromROSMsg(pcd, *map_points_ptr);


  pclndt_.setTransformationEpsilon(0.05);
  pclndt_.setStepSize(0.1);
  pclndt_.setResolution(2.0);
  pclndt_.setMaximumIterations(30);
  pclndt_.setInputTarget(map_points_ptr); 

#ifndef USE_IKFOM
  VD(DIM_STATE)
  solution;
  MD(DIM_STATE, DIM_STATE)
  G, H_T_H, I_STATE;
  V3D rot_add, t_add;
  StatesGroup state_propagat;
  PointType pointOri, pointSel, coeff;
#endif
  // PointCloudXYZI::Ptr corr_normvect(new PointCloudXYZI(100000, 1));
  int effect_feat_num = 0, frame_num = 0;
  double deltaT, deltaR, aver_time_consu = 0, aver_time_icp = 0,
                         aver_time_match = 0, aver_time_solve = 0,
                         aver_time_const_H_time = 0;

  FOV_DEG = (fov_deg + 10.0) > 179.9 ? 179.9 : (fov_deg + 10.0);
  HALF_FOV_COS = cos((FOV_DEG) * 0.5 * PI_M / 180.0);

  downSizeFilterSurf.setLeafSize(filter_size_surf_min, filter_size_surf_min,
                                 filter_size_surf_min);



  shared_ptr<ImuProcess> p_imu(new ImuProcess());
  
  extT << VEC_FROM_ARRAY(extrinT);
  extR << MAT_FROM_ARRAY(extrinR);
  gilever << VEC_FROM_ARRAY(ginslever);


  p_imu->set_extrinsic(extT, extR);
  p_imu->set_gyr_cov_scale(V3D(gyr_cov, gyr_cov, gyr_cov));
  p_imu->set_acc_cov_scale(V3D(acc_cov, acc_cov, acc_cov));
  p_imu->set_gyr_bias_cov(V3D(0.0001, 0.0001, 0.0001));
  p_imu->set_acc_bias_cov(V3D(0.0001, 0.0001, 0.0001));
  p_imu->set_imu_init_frame_num(imu_int_frame);

  if (!imu_en)
    p_imu->disable_imu();

#ifndef USE_IKFOM
  G.setZero();
  H_T_H.setZero();
  I_STATE.setIdentity();
#endif


  /*** debug record ***/
  FILE *fp;
  string pos_log_dir = root_dir + "/Log/pos_log.txt";
  fp = fopen(pos_log_dir.c_str(), "w");
  fout_img_pos.open(string(string(ROOT_DIR) + "PCD/img_pos.json"), ios::out);
  fout_pcd_pos.open(string(string(ROOT_DIR) + "PCD/scans_pos.json"), ios::out);
  fout_pre.open(DEBUG_FILE_DIR("mat_pre.txt"), ios::out);
  fout_out.open(DEBUG_FILE_DIR("mat_out.txt"), ios::out);
  fout_dbg.open(DEBUG_FILE_DIR("dbg.txt"), ios::out);

  if (fout_pre && fout_out)
    cout << "~~~~" << ROOT_DIR << " file opened" << endl;
  else
    cout << "~~~~" << ROOT_DIR << " doesn't exist" << endl;



#ifdef USE_ikdforest
  ikdforest.Set_balance_criterion_param(0.6);
  ikdforest.Set_delete_criterion_param(0.5);
#endif


  signal(SIGINT, SigHandle);
  ros::Rate rate(5000);
  bool status = ros::ok();

  while (status)
  {
    if (flg_exit)
      break;
    ros::spinOnce();
    if (!sync_packages(LidarMeasures))
    {
      status = ros::ok();
      // cv::waitKey(1);
      rate.sleep();
      continue;
    }

    /*** Packaged got ***/
    if (flg_reset)
    {
      ROS_WARN("reset when rosbag play back");
      p_imu->Reset();
      flg_reset = false;
      continue;
    }

    if (!is_first_frame)
    {
      first_lidar_time = LidarMeasures.lidar_beg_time;
      p_imu->first_lidar_time = first_lidar_time;
      is_first_frame = true;
      cout << "FIRST LIDAR FRAME!" << endl;
    }

    double t0, t1, t2, t3, t4, t5, match_start, solve_start, svd_time;

    match_time = kdtree_search_time = kdtree_search_counter = solve_time =
        solve_const_H_time = svd_time = 0;
    t0 = omp_get_wtime();

    imu_vector.clear();
    p_imu->Process2(LidarMeasures, state, feats_undistort, imu_vector);
    state_propagat = state;
 
    if (!init_map_wx)
    {
      if (!gtQueue.empty())
      {
        odometryMsgToAffine3f(gtQueue.front(), affine_imu_to_map); // 转到了IMU坐标系
        gtQueue.pop_front();
        Eigen::Affine3d temp_init_m = affine_imu_to_map.cast<double>();
        state.rot_end = temp_init_m.rotation();
        state.pos_end = temp_init_m.translation();
        cout << "init state pos:\n"
             << state.pos_end << endl;
        init_map_wx = true;

      }
      else
      {
        continue;
      }

    }
    else
    {
      initial_pose_matrix = pre_trans * delta_trans;
    }

    Eigen::Matrix3d rot_varpre = state.cov.block<3, 3>(0, 0);
    Eigen::Matrix3d t_varpre = state.cov.block<3, 3>(3, 3);
    double rot1 = RAD2DEG(sqrt(rot_varpre(0, 0)));
    double rot2 = RAD2DEG(sqrt(rot_varpre(1, 1)));
    double rot3 = RAD2DEG(sqrt(rot_varpre(2, 2)));
    double trans1 = (sqrt(t_varpre(0, 0)));
    double trans2 = (sqrt(t_varpre(1, 1)));
    double trans3 = (sqrt(t_varpre(2, 2)));
    std::cout << std::setprecision(5) << "pre rot sigma:" << rot1 << " " << rot2 << " " << rot3 << " "
              << ",t sigma:" << trans1 << " " << trans2 << " " << trans3 << std::endl;

    euler_cur = RotMtoEuler(state.rot_end);

    fout_pre << setw(20) << LidarMeasures.lidar_beg_time - first_lidar_time
             << " " << euler_cur.transpose() * 180 / M_PI << " "
             << state.pos_end.transpose() << " " << state.vel_end.transpose()
             << " " << state.bias_g.transpose() << " "
             << state.bias_a.transpose() << " " << state.gravity.transpose()
             << endl;
    if (feats_undistort->empty() || (feats_undistort == nullptr))
    {
      cout << " No point!!!" << endl;
      continue;
    }

    flg_EKF_inited = (LidarMeasures.lidar_beg_time - first_lidar_time) < INIT_TIME ? false : true;


    feats_down_body = feats_undistort;
    int undistor_size = feats_undistort->size();

    pcl::PointCloud<pcl::PointXYZI> ori_cloud;
    for (size_t i = 0; i < feats_undistort->size(); i++)
    {
      pcl::PointXYZI pi;
      pi.x = feats_undistort->points[i].x;
      pi.y = feats_undistort->points[i].y;
      pi.z = feats_undistort->points[i].z;
      pi.intensity = feats_undistort->points[i].intensity;
      ori_cloud.push_back(pi);
    }
    down_sampling_voxel(ori_cloud, filter_size_surf_min);
    feats_down_body->clear();
    pcl::PointCloud<pcl::PointXYZ>::Ptr ori_cloud_ndt(new pcl::PointCloud<pcl::PointXYZ>);

    for (size_t i = 0; i < ori_cloud.size(); i++)
    {
      pcl::PointXYZINormal pi;
      pi.x = ori_cloud.points[i].x;
      pi.y = ori_cloud.points[i].y;
      pi.z = ori_cloud.points[i].z;
      pi.intensity = ori_cloud.points[i].intensity;
      feats_down_body->push_back(pi);

      pcl::PointXYZ pi_2;
      V3D point_this(ori_cloud.points[i].x, ori_cloud.points[i].y,
                     ori_cloud.points[i].z);
      point_this = extR * point_this + extT;
      pi_2.x = point_this[0];
      pi_2.y = point_this[1];
      pi_2.z = point_this[2];
      ori_cloud_ndt->push_back(pi_2);
    }

    std::cout << "!!!!!!!!! voxel grid filter:" << undistor_size << " -> " << feats_down_body->size() << std::endl;
    feats_down_size = feats_down_body->points.size();


    /*** ICP and iterated Kalman filter update ***/
    normvec->resize(feats_down_size);
    feats_down_world->clear();
    feats_down_world->resize(feats_down_size);
    // vector<double> res_last(feats_down_size, 1000.0); // initial //
    res_last.resize(feats_down_size, 1000.0);
    point_selected_surf.resize(feats_down_size, true);
    pointSearchInd_surf.resize(feats_down_size);
    Nearest_Points.resize(feats_down_size);

    StatesGroup state_gnss = state;

    Eigen::Matrix3d rot_varlidar = state.cov.block<3, 3>(0, 0);
    Eigen::Matrix3d t_varlidar = state.cov.block<3, 3>(3, 3);
    double rot1l = RAD2DEG(sqrt(rot_varlidar(0, 0)));
    double rot2l = RAD2DEG(sqrt(rot_varlidar(1, 1)));
    double rot3l = RAD2DEG(sqrt(rot_varlidar(2, 2)));
    double trans1l = (sqrt(t_varlidar(0, 0)));
    double trans2l = (sqrt(t_varlidar(1, 1)));
    double trans3l = (sqrt(t_varlidar(2, 2)));
    std::cout << std::setprecision(5) << "pre lidar rot sigma:" << rot1l << " " << rot2l << " " << rot3l << " "
              << ",lidar t sigma:" << trans1l << " " << trans2l << " " << trans3l << std::endl;

    Eigen::Matrix4d predict_state;
    predict_state.setZero();
    predict_state.block(0, 0, 3, 3) = state.rot_end;
    predict_state.block(0, 3, 3, 1) = state.pos_end;
    predict_state(3, 3) = 1;

  // pclndt
    pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    Sophus::SE3d guess_imu = Sophus::SE3d::fitToSE3(predict_state);
    Eigen::Matrix4f result_pose_matrix = guess_imu.matrix().cast<float>(); // 将SE3对象guess转换为Matrix4f格式
    pclndt_.setInputSource(ori_cloud_ndt);
    pclndt_.align(*output_cloud, result_pose_matrix);
    double fitness_score;
    fitness_score = pclndt_.getFitnessScore(); // 两片点云之间的误差 
    result_pose_matrix = pclndt_.getFinalTransformation();
    Eigen::Matrix4d map_op = result_pose_matrix.matrix().cast<double>();
    V3D observe_pose = map_op.block(0, 3, 3, 1);
    M3D observe_rot = map_op.block(0, 0, 3, 3);

    cout << "observe_pose:\n"
         << observe_pose << endl;


    V3D euler_cur_map = RotMtoEuler(observe_rot);


    double total_distance_ndt = (observe_pose - position_last_ndt).norm();

    printf(REDPURPLE "[ndt move  distance] ");
    cout << " total_distance_ndt:  " << total_distance_ndt << endl;
    cout << RESET;

    loose_observe(state, observe_pose, observe_rot, 0.1, 0.05); 

    Eigen::Matrix4d Result_state;
    Result_state.setZero();
    Result_state.block(0, 0, 3, 3) = state.rot_end;
    Result_state.block(0, 3, 3, 1) = state.pos_end;
    Result_state(3, 3) = 1;
    Sophus::SE3d result_kf = Sophus::SE3d::fitToSE3(Result_state);
    Eigen::Matrix4f result_kf2 = result_kf.matrix().cast<float>(); 
    delta_trans = pre_trans.inverse() * result_kf2;
    pre_trans = result_kf2;

    pcl::PointCloud<pcl::PointXYZ>::Ptr sensor_points_mapTF_ptr(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::transformPointCloud(
        *ori_cloud_ndt, *sensor_points_mapTF_ptr, result_kf2);
    sensor_msgs::PointCloud2 sensor_points_mapTF_msg;
    pcl::toROSMsg(*sensor_points_mapTF_ptr, sensor_points_mapTF_msg);
    sensor_points_mapTF_msg.header.stamp = ros::Time::now();
    sensor_points_mapTF_msg.header.frame_id = "camera_init";
    sensor_aligned_pose_pub_.publish(sensor_points_mapTF_msg);

    position_last_ndt = state.pos_end;

    // double t_update_end = omp_get_wtime();
    /******* Publish odometry *******/
    euler_cur = RotMtoEuler(state.rot_end);
    cout << "euler_cur:\n"
         << euler_cur << endl;
    geoQuat = tf::createQuaternionMsgFromRollPitchYaw(euler_cur(0), euler_cur(1), euler_cur(2));
    publish_odometry(pubOdomAftMapped);

    Eigen::Matrix3d rot_var = state.cov.block<3, 3>(0, 0);
    Eigen::Matrix3d t_var = state.cov.block<3, 3>(3, 3);
    double r1_e = RAD2DEG(sqrt(rot_var(0, 0)));
    double r2_e = RAD2DEG(sqrt(rot_var(1, 1)));
    double r3_e = RAD2DEG(sqrt(rot_var(2, 2)));
    double t1_e = (sqrt(t_var(0, 0)));
    double t2_e = (sqrt(t_var(1, 1)));
    double t3_e = (sqrt(t_var(2, 2)));
    std::cout << std::setprecision(5) << "out rot sigma:" << r1_e << " " << r2_e << " " << r3_e << " "
              << ",t sigma:" << t1_e << " " << t2_e << " " << t3_e << std::endl;


    kdtree_incremental_time = t5 - t3 + readd_time;

    int size_un = feats_undistort->points.size();
    PointCloudXYZI::Ptr laserCloudWorld_merge(new PointCloudXYZI(size_un, 1));
    for (int i = 0; i < size_un; i++)
    {
      RGBpointBodyToWorld(&feats_down_body->points[i],
                          &laserCloudWorld_merge->points[i]);
    }

    *pcl_wait_pub = *laserCloudWorld_merge;

    publish_path(pubPath);
      euler_cur = RotMtoEuler(state.rot_end);

      fout_out << std::setw(20) << LidarMeasures.lidar_beg_time - first_lidar_time
               << " " << euler_cur.transpose() * 180 / M_PI << " "
               << state.pos_end.transpose() << " " << state.vel_end.transpose()
               << " " << state.bias_g.transpose() << " "
               << state.bias_a.transpose() << " " << state.gravity.transpose()
               << " " << feats_undistort->points.size() << endl;

    double t_end_all = omp_get_wtime();


  }

  return 0;
}