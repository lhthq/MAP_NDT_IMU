#ifndef IKFOM_TYPES_H
#define IKFOM_TYPES_H

#include <IKFoM_toolkit/esekfom/esekfom.hpp>

/*** For MTK ***/
typedef MTK::vect<3, double> vect3;
typedef MTK::SO3<double> SO3_data;
typedef MTK::S2<double, 98090, 10000, 1> S2; 
typedef MTK::vect<1, double> vect1;
typedef MTK::vect<2, double> vect2;

MTK_BUILD_MANIFOLD(state_ikfom,
((vect3, pos))
((SO3_data, rot))
((SO3_data, offset_R_L_I))
((vect3, offset_T_L_I))
((vect3, vel))
((vect3, bg))
((vect3, ba))
((S2, grav))
);

MTK_BUILD_MANIFOLD(input_ikfom,
((vect3, acc))
((vect3, gyro))
);

MTK_BUILD_MANIFOLD(process_noise_ikfom,
((vect3, ng))
((vect3, na))
((vect3, nbg))
((vect3, nba))
);


#endif