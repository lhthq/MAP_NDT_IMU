#include "IMU_Processing.h"

const bool time_list(PointType &x, PointType &y)
{
  return (x.curvature < y.curvature);
}

ImuProcess::ImuProcess()
    : b_first_frame_(true), imu_need_init_(true), start_timestamp_(-1)
{
  init_iter_num = 1;


  for (size_t i = 0; i < 3; i++)
  {
    neh[i] = 0.0;
    cond[i] = 0.0;
    init_pose[i] = 0.0;
  }

  for (int ii = 0; ii < 4; ii++)
  {
    chose_anchor[ii] = 0;
  }

  cov_acc = V3D(0.1, 0.1, 0.1);
  cov_gyr = V3D(0.1, 0.1, 0.1);

  cov_bias_gyr = V3D(0.1, 0.1, 0.1);
  cov_bias_acc = V3D(0.1, 0.1, 0.1);
  mean_acc = V3D(0, 0, -1.0);
  mean_gyr = V3D(0, 0, 0);
  angvel_last = Zero3d;
  acc_s_last = Zero3d;
  Lid_offset_to_IMU = Zero3d;
  Lid_rot_to_IMU = Eye3d;
  last_imu_.reset(new sensor_msgs::Imu());
}

ImuProcess::~ImuProcess() {}

void ImuProcess::Reset()
{
  ROS_WARN("Reset ImuProcess");
  mean_acc = V3D(0, 0, -1.0);
  mean_gyr = V3D(0, 0, 0);
  angvel_last = Zero3d;
  imu_need_init_ = true;
  start_timestamp_ = -1;
  init_iter_num = 1;
  v_imu_.clear();
  IMUpose.clear();
  last_imu_.reset(new sensor_msgs::Imu());
  cur_pcl_un_.reset(new PointCloudXYZI());
}

void ImuProcess::disable_imu()
{
  cout << "IMU disabled !!!!!" << endl;
  imu_en = false;
  imu_need_init_ = false;
}

void ImuProcess::push_update_state(double offs_t, StatesGroup state)
{

  V3D acc_tmp = acc_s_last, angvel_tmp = angvel_last, vel_imu(state.vel_end), pos_imu(state.pos_end);
  M3D R_imu(state.rot_end);
  IMUpose.push_back(set_pose6d(offs_t, acc_tmp, angvel_tmp, vel_imu, pos_imu, R_imu));
}

void ImuProcess::set_extrinsic(const MD(4, 4) & T)
{
  Lid_offset_to_IMU = T.block<3, 1>(0, 3);
  Lid_rot_to_IMU = T.block<3, 3>(0, 0);
}

void ImuProcess::set_extrinsic(const V3D &transl)
{
  Lid_offset_to_IMU = transl;
  Lid_rot_to_IMU.setIdentity();
}

void ImuProcess::set_extrinsic(const V3D &transl, const M3D &rot)
{
  Lid_offset_to_IMU = transl;
  Lid_rot_to_IMU = rot;
}

void ImuProcess::set_gyr_cov_scale(const V3D &scaler)
{
  cov_gyr = scaler;
}

void ImuProcess::set_acc_cov_scale(const V3D &scaler)
{
  cov_acc = scaler;
}

void ImuProcess::set_gyr_bias_cov(const V3D &b_g)
{
  cov_bias_gyr = b_g;
}

void ImuProcess::set_acc_bias_cov(const V3D &b_a)
{
  cov_bias_acc = b_a;
}

void ImuProcess::set_imu_init_frame_num(const int &num)
{
  MAX_INI_COUNT = num;
}


void ImuProcess::IMU_init(const MeasureGroup &meas, StatesGroup &state_inout, int &N)
{

  ROS_INFO("IMU Initializing: %.1f %%", double(N) / MAX_INI_COUNT * 100);
  V3D cur_acc, cur_gyr;
  V3D fir_iter;
  bool att_result = false;
  double yaw = -300;

  if (b_first_frame_)
  {
    Reset();
    N = 1;

    b_first_frame_ = false;

    const auto &imu_acc = meas.imu.front()->linear_acceleration;
    const auto &gyr_acc = meas.imu.front()->angular_velocity;
    mean_acc << imu_acc.x, imu_acc.y, imu_acc.z;
    mean_gyr << gyr_acc.x, gyr_acc.y, gyr_acc.z;
  }


  for (const auto &imu : meas.imu)
  {
    const auto &imu_acc = imu->linear_acceleration;
    const auto &gyr_acc = imu->angular_velocity;
    cur_acc << imu_acc.x, imu_acc.y, imu_acc.z;
    cur_gyr << gyr_acc.x, gyr_acc.y, gyr_acc.z;

    mean_acc += (cur_acc - mean_acc) / N;
    mean_gyr += (cur_gyr - mean_gyr) / N;

    N++;
  }


  state_inout.gravity = -mean_acc / mean_acc.norm() * G_m_s2;
  // state_inout.rot_end = Eye3d; // Exp(mean_acc.cross(V3D(0, 0, -1 / scale_gravity)));
  state_inout.bias_g = mean_gyr; // Zero3d; // mean_gyr;
  last_imu_ = meas.imu.back();
}

void ImuProcess::Forward(const MeasureGroup &meas, StatesGroup &state_inout, double pcl_beg_time, double end_time)
{
  /*** add the imu of the last frame-tail to the of current frame-head ***/
  auto v_imu = meas.imu;
  v_imu.push_front(last_imu_);
  const double &imu_beg_time = v_imu.front()->header.stamp.toSec();
  const double &imu_end_time = v_imu.back()->header.stamp.toSec();


  if (IMUpose.empty())
  {
    IMUpose.push_back(set_pose6d(0.0, acc_s_last, angvel_last, state_inout.vel_end, state_inout.pos_end, state_inout.rot_end));
  }

  /*** forward propagation at each imu point ***/
  V3D acc_imu = acc_s_last, angvel_avr = angvel_last, acc_avr, vel_imu(state_inout.vel_end), pos_imu(state_inout.pos_end);
  M3D R_imu(state_inout.rot_end);
  //  last_state = state_inout;
  MD(DIM_STATE, DIM_STATE)
  F_x, cov_w;

  double dt = 0;
  for (auto it_imu = v_imu.begin(); it_imu < (v_imu.end() - 1); it_imu++)
  {
    auto &&head = *(it_imu);
    auto &&tail = *(it_imu + 1);

    if (tail->header.stamp.toSec() < last_lidar_end_time_)
      continue;

    angvel_avr << 0.5 * (head->angular_velocity.x + tail->angular_velocity.x),
        0.5 * (head->angular_velocity.y + tail->angular_velocity.y),
        0.5 * (head->angular_velocity.z + tail->angular_velocity.z);


    acc_avr << 0.5 * (head->linear_acceleration.x + tail->linear_acceleration.x),
        0.5 * (head->linear_acceleration.y + tail->linear_acceleration.y),
        0.5 * (head->linear_acceleration.z + tail->linear_acceleration.z);
    last_acc = acc_avr;
    last_ang = angvel_avr;
    // #ifdef DEBUG_PRINT
    fout_imu << setw(10) << head->header.stamp.toSec() - first_lidar_time << " " << angvel_avr.transpose() << " " << acc_avr.transpose() << endl;
    // #endif

    angvel_avr -= state_inout.bias_g;
    acc_avr = acc_avr * G_m_s2 / mean_acc.norm() - state_inout.bias_a;

    if (head->header.stamp.toSec() < last_lidar_end_time_)
    {
      dt = tail->header.stamp.toSec() - last_lidar_end_time_;
    }
    else
    {
      dt = tail->header.stamp.toSec() - head->header.stamp.toSec();
    }
    // cout<<setw(20)<<"dt: "<<dt<<endl;
    /* covariance propagation */
    M3D acc_avr_skew;
    M3D Exp_f = Exp(angvel_avr, dt);
    acc_avr_skew << SKEW_SYM_MATRX(acc_avr);

    F_x.setIdentity();
    cov_w.setZero();

    F_x.block<3, 3>(0, 0) = Exp(angvel_avr, -dt);
    F_x.block<3, 3>(0, 9) = -Eye3d * dt;
    // F_x.block<3,3>(3,0)  = R_imu * off_vel_skew * dt;
    F_x.block<3, 3>(3, 6) = Eye3d * dt;
    F_x.block<3, 3>(6, 0) = -R_imu * acc_avr_skew * dt;
    F_x.block<3, 3>(6, 12) = -R_imu * dt;
    F_x.block<3, 3>(6, 15) = Eye3d * dt;

    cov_w.block<3, 3>(0, 0).diagonal() = cov_gyr * dt * dt;
    cov_w.block<3, 3>(6, 6) = R_imu * cov_acc.asDiagonal() * R_imu.transpose() * dt * dt;
    cov_w.block<3, 3>(9, 9).diagonal() = cov_bias_gyr * dt * dt;   // bias gyro covariance
    cov_w.block<3, 3>(12, 12).diagonal() = cov_bias_acc * dt * dt; // bias acc covariance

    state_inout.cov = F_x * state_inout.cov * F_x.transpose() + cov_w;

    /* propogation of IMU attitude */
    R_imu = R_imu * Exp_f;

    /* Specific acceleration (global frame) of IMU */
    acc_imu = R_imu * acc_avr + state_inout.gravity;

    /* propogation of IMU */
    pos_imu = pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt;

    /* velocity of IMU */
    vel_imu = vel_imu + acc_imu * dt;

    /* save the poses at each IMU measurements */
    angvel_last = angvel_avr;
    acc_s_last = acc_imu;
    double &&offs_t = tail->header.stamp.toSec() - pcl_beg_time;
    IMUpose.push_back(set_pose6d(offs_t, acc_imu, angvel_avr, vel_imu, pos_imu, R_imu));
  }

  /*** calculated the pos and attitude prediction at the frame-end ***/
  double note = end_time > imu_end_time ? 1.0 : -1.0;
  dt = note * (end_time - imu_end_time);
  state_inout.vel_end = vel_imu + note * acc_imu * dt;
  state_inout.rot_end = R_imu * Exp(V3D(note * angvel_avr), dt);
  state_inout.pos_end = pos_imu + note * vel_imu * dt + note * 0.5 * acc_imu * dt * dt;

  last_imu_ = v_imu.back();
  last_lidar_end_time_ = end_time;

  // auto pos_liD_e = state_inout.pos_end + state_inout.rot_end * Lid_offset_to_IMU;
  // auto R_liD_e   = state_inout.rot_end * Lidar_R_to_IMU;


}

void ImuProcess::Forward_without_imu(LidarMeasureGroup &meas, StatesGroup &state_inout, PointCloudXYZI &pcl_out)
{
  const double &pcl_beg_time = meas.lidar_beg_time;

  /*** sort point clouds by offset time ***/
  pcl_out = *(meas.lidar);
  // sort(pcl_out->points.begin(), pcl_out->points.end(), time_list);
  const double &pcl_end_time =
      pcl_beg_time + pcl_out.points.back().curvature / double(1000);
  // V3D acc_imu, angvel_avr, acc_avr, vel_imu(state_inout.vel_end),
  //     pos_imu(state_inout.pos_end);
  // M3D R_imu(state_inout.rot_end);
  meas.last_update_time = pcl_end_time;
  MD(DIM_STATE, DIM_STATE)
  F_x, cov_w;
  double dt = 0;

  if (b_first_frame_)
  {
    dt = 0.1;
    b_first_frame_ = false;
  }
  else
  {
    dt = pcl_beg_time - time_last_scan;
  }

  time_last_scan = pcl_beg_time;

  M3D Exp_f = Exp(state_inout.bias_g, dt);

  F_x.setIdentity();
  cov_w.setZero();

  F_x.block<3, 3>(0, 0) = Exp(state_inout.bias_g, -dt);
  F_x.block<3, 3>(0, 9) = Eye3d * dt;
  F_x.block<3, 3>(3, 6) = Eye3d * dt;
  // F_x.block<3, 3>(6, 0)  = - R_imu * acc_avr_skew * dt;
  // F_x.block<3, 3>(6, 12) = - R_imu * dt;
  // F_x.block<3, 3>(6, 15) = Eye3d * dt;

  cov_w.block<3, 3>(9, 9).diagonal() = cov_gyr * dt * dt; // for omega in constant model
  cov_w.block<3, 3>(6, 6).diagonal() = cov_acc * dt * dt; // for velocity in constant model


  // std::cout << "before propagete:" << state_inout.cov.diagonal().transpose()
  //           << std::endl;
  state_inout.cov = F_x * state_inout.cov * F_x.transpose() + cov_w;
  // std::cout << "cov_w:" << cov_w.diagonal().transpose() << std::endl;
  // std::cout << "after propagete:" << state_inout.cov.diagonal().transpose()
  //           << std::endl;
  state_inout.rot_end = state_inout.rot_end * Exp_f;
  state_inout.pos_end = state_inout.pos_end + state_inout.vel_end * dt;
}

void ImuProcess::Backward(const LidarMeasureGroup &lidar_meas, StatesGroup &state_inout, PointCloudXYZI &pcl_out)
{
  /*** undistort each lidar point (backward propagation) ***/
  M3D R_imu;
  V3D acc_imu, angvel_avr, vel_imu, pos_imu;
  double dt;
  auto pos_liD_e = state_inout.pos_end + state_inout.rot_end * Lid_offset_to_IMU;
  auto it_pcl = pcl_out.points.end() - 1;
  for (auto it_kp = IMUpose.end() - 1; it_kp != IMUpose.begin(); it_kp--)
  {
    auto head = it_kp - 1;
    auto tail = it_kp;
    R_imu << MAT_FROM_ARRAY(head->rot);
    acc_imu << VEC_FROM_ARRAY(head->acc);
    // cout<<"head imu acc: "<<acc_imu.transpose()<<endl;
    vel_imu << VEC_FROM_ARRAY(head->vel);
    pos_imu << VEC_FROM_ARRAY(head->pos);
    angvel_avr << VEC_FROM_ARRAY(head->gyr);
    for (; it_pcl->curvature / double(1000) > head->offset_time; it_pcl--)
    {
      dt = it_pcl->curvature / double(1000) - head->offset_time;

      /* Transform to the 'end' frame, using only the rotation
       * Note: Compensation direction is INVERSE of Frame's moving direction
       * So if we want to compensate a point at timestamp-i to the frame-e
       * P_compensate = R_imu_e ^ T * (R_i * P_i + T_ei) where T_ei is represented in global frame */
      M3D R_i(R_imu * Exp(angvel_avr, dt));
      V3D T_ei(pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt + R_i * Lid_offset_to_IMU - pos_liD_e);

      V3D P_i(it_pcl->x, it_pcl->y, it_pcl->z);
      V3D P_compensate = state_inout.rot_end.transpose() * (R_i * P_i + T_ei);

      /// save Undistorted points and their rotation
      it_pcl->x = P_compensate(0);
      it_pcl->y = P_compensate(1);
      it_pcl->z = P_compensate(2);

      if (it_pcl == pcl_out.points.begin())
        break;
    }
  }
}

void ImuProcess::Process(LidarMeasureGroup &lidar_meas, StatesGroup &stat, PointCloudXYZI::Ptr cur_pcl_un_)
{
  double t1, t2, t3;
  // t1 = omp_get_wtime();
  ROS_ASSERT(lidar_meas.lidar != nullptr);
  MeasureGroup meas = lidar_meas.measures.back();

  if (imu_need_init_)
  {
    if (meas.imu.empty())
    {
      return;
    };
    /// The very first lidar frame
    IMU_init(meas, stat, init_iter_num);

    imu_need_init_ = true;

    last_imu_ = meas.imu.back();

    if (init_iter_num > MAX_INI_COUNT)
    {
      cov_acc *= pow(G_m_s2 / mean_acc.norm(), 2);
      imu_need_init_ = false;
      ROS_INFO("IMU Initials: Gravity: %.4f %.4f %.4f %.4f; acc covarience: %.8f %.8f %.8f; gry covarience: %.8f %.8f %.8f",
               stat.gravity[0], stat.gravity[1], stat.gravity[2], mean_acc.norm(), cov_acc[0], cov_acc[1], cov_acc[2], cov_gyr[0], cov_gyr[1], cov_gyr[2]);

      // cout<<"mean acc: "<<mean_acc<<" acc measures in word frame:"<<state.rot_end.transpose()*mean_acc<<endl;
      fout_imu.open(DEBUG_FILE_DIR("imu.txt"), ios::out);
    }

    return;
  }

  /// Undistort points： the first point is assummed as the base frame
  /// Compensate lidar points with IMU rotation (with only rotation now)
  if (lidar_meas.is_lidar_end)
  {
    /*** sort point clouds by offset time ***/
    *cur_pcl_un_ = *(lidar_meas.lidar);
    sort(cur_pcl_un_->points.begin(), cur_pcl_un_->points.end(), time_list);
    const double &pcl_beg_time = lidar_meas.lidar_beg_time;
    const double &pcl_end_time = pcl_beg_time + lidar_meas.lidar->points.back().curvature / double(1000);
    if (imu_en)
    {
      Forward(meas, stat, pcl_beg_time, pcl_end_time);
      Backward(lidar_meas, stat, *cur_pcl_un_);
      last_lidar_end_time_ = pcl_end_time;
      IMUpose.clear();
    }
    else
    {
      Forward_without_imu(lidar_meas, stat, *cur_pcl_un_);
    }

  }
  else
  {
    const double &pcl_beg_time = lidar_meas.lidar_beg_time;
    const double &img_end_time = pcl_beg_time + meas.img_offset_time;
    Forward(meas, stat, pcl_beg_time, img_end_time);
  }

}

void ImuProcess::UndistortPcl(LidarMeasureGroup &lidar_meas, StatesGroup &state_inout, PointCloudXYZI &pcl_out, vector<Pose6D> &imu_vector)
{
  static int imu_count = 0;
  /*** add the imu of the last frame-tail to the of current frame-head ***/
  MeasureGroup meas;
  meas = lidar_meas.measures.back();
  // cout<<"meas.imu.size: "<<meas.imu.size()<<endl;
  auto v_imu = meas.imu;
  v_imu.push_front(last_imu_);
  const double &imu_beg_time = v_imu.front()->header.stamp.toSec();
  const double &imu_end_time = v_imu.back()->header.stamp.toSec();
  const double pcl_beg_time = MAX(lidar_meas.lidar_beg_time, lidar_meas.last_update_time);
  // const double &pcl_beg_time = meas.lidar_beg_time;

  /*** sort point clouds by offset time ***/
  pcl_out.clear();
  
  const double pcl_end_time = lidar_meas.lidar_beg_time + lidar_meas.lidar->points.back().curvature / double(1000); // lidar_meas.is_lidar_end ? lidar_meas.lidar_beg_time + lidar_meas.lidar->points.back().curvature / double(1000) : lidar_meas.lidar_beg_time + lidar_meas.measures.back().img_offset_time;

  pcl_out = *(lidar_meas.lidar);

  /*** Initialize IMU pose ***/
  IMUpose.clear();
  // IMUpose.push_back(set_pose6d(0.0, Zero3d, Zero3d, state.vel_end, state.pos_end, state.rot_end));
  IMUpose.push_back(set_pose6d(0.0, acc_s_last, angvel_last, state_inout.vel_end, state_inout.pos_end, state_inout.rot_end));

  /*** forward propagation at each imu point ***/
  V3D acc_imu(acc_s_last), angvel_avr(angvel_last), acc_avr, vel_imu(state_inout.vel_end), pos_imu(state_inout.pos_end);
  M3D R_imu(state_inout.rot_end);
  MD(DIM_STATE, DIM_STATE)
  F_x, cov_w;

  double dt = 0;
  for (auto it_imu = v_imu.begin(); it_imu != v_imu.end() - 1; it_imu++)
  {

    imu_count++;

    auto &&head = *(it_imu);
    auto &&tail = *(it_imu + 1);

    if (tail->header.stamp.toSec() < last_lidar_end_time_)
      continue;

    angvel_avr << 0.5 * (head->angular_velocity.x + tail->angular_velocity.x),
        0.5 * (head->angular_velocity.y + tail->angular_velocity.y),
        0.5 * (head->angular_velocity.z + tail->angular_velocity.z);

    // angvel_avr<<tail->angular_velocity.x, tail->angular_velocity.y, tail->angular_velocity.z;

    acc_avr << 0.5 * (head->linear_acceleration.x + tail->linear_acceleration.x),
        0.5 * (head->linear_acceleration.y + tail->linear_acceleration.y),
        0.5 * (head->linear_acceleration.z + tail->linear_acceleration.z);

    // #ifdef DEBUG_PRINT
    fout_imu << setw(10) << head->header.stamp.toSec() - first_lidar_time << " " << angvel_avr.transpose() << " " << acc_avr.transpose() << endl;
    // #endif

    angvel_avr -= state_inout.bias_g;
    acc_avr = acc_avr * G_m_s2 / mean_acc.norm() - state_inout.bias_a;

    if (head->header.stamp.toSec() < last_lidar_end_time_)
    {
      dt = tail->header.stamp.toSec() - last_lidar_end_time_;
    }
    else
    {
      dt = tail->header.stamp.toSec() - head->header.stamp.toSec();
    }

    /* covariance propagation */
    M3D acc_avr_skew;
    M3D Exp_f = Exp(angvel_avr, dt);
    acc_avr_skew << SKEW_SYM_MATRX(acc_avr);

    F_x.setIdentity();
    cov_w.setZero();

    F_x.block<3, 3>(0, 0) = Exp(angvel_avr, -dt);
    F_x.block<3, 3>(0, 9) = -Eye3d * dt; //
    // F_x.block<3,3>(3,0)  = R_imu * off_vel_skew * dt;
    F_x.block<3, 3>(3, 6) = Eye3d * dt;
    F_x.block<3, 3>(6, 0) = -R_imu * acc_avr_skew * dt;
    F_x.block<3, 3>(6, 12) = -R_imu * dt;
    F_x.block<3, 3>(6, 15) = Eye3d * dt;


    cov_w.block<3, 3>(0, 0).diagonal() = cov_gyr * dt * dt;
    cov_w.block<3, 3>(6, 6) = R_imu * cov_acc.asDiagonal() * R_imu.transpose() * dt * dt;
    cov_w.block<3, 3>(9, 9).diagonal() = cov_bias_gyr * dt * dt;   // bias gyro covariance
    cov_w.block<3, 3>(12, 12).diagonal() = cov_bias_acc * dt * dt; // bias acc covariance



    state_inout.cov = F_x * state_inout.cov * F_x.transpose() + cov_w;


    /* propogation of IMU attitude */
    R_imu = R_imu * Exp_f;

    /* Specific acceleration (global frame) of IMU */
    acc_imu = R_imu * acc_avr + state_inout.gravity;

    /* propogation of IMU */
    pos_imu = pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt;

    /* velocity of IMU */
    vel_imu = vel_imu + acc_imu * dt;

    /* save the poses at each IMU measurements */
    angvel_last = angvel_avr;
    acc_s_last = acc_imu;
    double &&offs_t = tail->header.stamp.toSec() - pcl_beg_time;
    // cout<<setw(20)<<"offset_t: "<<offs_t<<"tail->header.stamp.toSec(): "<<tail->header.stamp.toSec()<<endl;
    IMUpose.push_back(set_pose6d(offs_t, acc_imu, angvel_avr, vel_imu, pos_imu, R_imu));
    imu_vector.push_back(set_pose6d(offs_t, acc_imu, angvel_avr, vel_imu, pos_imu, R_imu));
  }

  /*** calculated the pos and attitude prediction at the frame-end ***/
  if (imu_end_time > pcl_beg_time)
  {
    double note = pcl_end_time > imu_end_time ? 1.0 : -1.0;
    dt = note * (pcl_end_time - imu_end_time);
    state_inout.vel_end = vel_imu + note * acc_imu * dt;
    state_inout.rot_end = R_imu * Exp(V3D(note * angvel_avr), dt);
    state_inout.pos_end = pos_imu + note * vel_imu * dt + note * 0.5 * acc_imu * dt * dt;
  }
  else
  {
    double note = pcl_end_time > pcl_beg_time ? 1.0 : -1.0;
    dt = note * (pcl_end_time - pcl_beg_time);
    state_inout.vel_end = vel_imu + note * acc_imu * dt;
    state_inout.rot_end = R_imu * Exp(V3D(note * angvel_avr), dt);
    state_inout.pos_end = pos_imu + note * vel_imu * dt + note * 0.5 * acc_imu * dt * dt;
  }

  last_imu_ = v_imu.back();
  last_lidar_end_time_ = pcl_end_time;

  if (pcl_out.points.size() < 1)
    return;
  /*** undistort each lidar point (backward propagation) ***/
  auto it_pcl = pcl_out.points.end() - 1;
  for (auto it_kp = IMUpose.end() - 1; it_kp != IMUpose.begin(); it_kp--)
  {
    auto head = it_kp - 1;
    auto tail = it_kp;
    R_imu << MAT_FROM_ARRAY(head->rot);
    acc_imu << VEC_FROM_ARRAY(head->acc);
    // cout<<"head imu acc: "<<acc_imu.transpose()<<endl;
    vel_imu << VEC_FROM_ARRAY(head->vel);
    pos_imu << VEC_FROM_ARRAY(head->pos);
    angvel_avr << VEC_FROM_ARRAY(head->gyr);

    for (; it_pcl->curvature / double(1000) > head->offset_time; it_pcl--)
    {
      dt = it_pcl->curvature / double(1000) - head->offset_time;

      /* Transform to the 'end' frame */
      M3D R_i(R_imu * Exp(angvel_avr, dt));
      V3D T_ei(pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt - state_inout.pos_end);

      V3D P_i(it_pcl->x, it_pcl->y, it_pcl->z);
      V3D P_compensate = Lid_rot_to_IMU.transpose() * (state_inout.rot_end.transpose() * (R_i * (Lid_rot_to_IMU * P_i + Lid_offset_to_IMU) + T_ei) - Lid_offset_to_IMU);

      /// save Undistorted points and their rotation
      it_pcl->x = P_compensate(0);
      it_pcl->y = P_compensate(1);
      it_pcl->z = P_compensate(2);

      if (it_pcl == pcl_out.points.begin())
        break;
    }
  }
  // cout << "[ IMU Process ]: undistort size: " << pcl_out.points.size() << endl;
}

void ImuProcess::Process2(LidarMeasureGroup &lidar_meas, StatesGroup &stat, PointCloudXYZI::Ptr cur_pcl_un_, vector<Pose6D> &imu_vector)
{
  double t1, t2, t3;
  // t1 = omp_get_wtime();
  ROS_ASSERT(lidar_meas.lidar != nullptr);
  if (!imu_en)
  {
    Forward_without_imu(lidar_meas, stat, *cur_pcl_un_);
    return;
  }

  MeasureGroup meas = lidar_meas.measures.back();

  if (imu_need_init_)
  {
    double pcl_end_time = lidar_meas.is_lidar_end ? lidar_meas.lidar_beg_time + lidar_meas.lidar->points.back().curvature / double(1000) : lidar_meas.lidar_beg_time + lidar_meas.measures.back().img_offset_time;
    lidar_meas.last_update_time = pcl_end_time;

    if (meas.imu.empty())
    {
      return;
    };
    /// The very first lidar frame      UWB->aid IMU initial    liuhong20230416
    IMU_init(meas, stat, init_iter_num);

    imu_need_init_ = true;

    last_imu_ = meas.imu.back();

    if (init_iter_num > MAX_INI_COUNT)
    {
      cov_acc *= pow(G_m_s2 / mean_acc.norm(), 2);
      imu_need_init_ = false;
      ROS_INFO("IMU Initials: Gravity: %.4f %.4f %.4f %.4f; acc covarience: %.8f %.8f %.8f; gry covarience: %.8f %.8f %.8f \n",
               stat.gravity[0], stat.gravity[1], stat.gravity[2], mean_acc.norm(), cov_acc[0], cov_acc[1], cov_acc[2], cov_gyr[0], cov_gyr[1], cov_gyr[2]);
      ROS_INFO("IMU Initials: ba covarience: %.8f %.8f %.8f; bg covarience: %.8f %.8f %.8f",
               cov_bias_acc[0], cov_bias_acc[1], cov_bias_acc[2], cov_bias_gyr[0], cov_bias_gyr[1], cov_bias_gyr[2]);
      fout_imu.open(DEBUG_FILE_DIR("imu.txt"), ios::out);
    }

    return;
  }

  UndistortPcl(lidar_meas, stat, *cur_pcl_un_, imu_vector);
}

