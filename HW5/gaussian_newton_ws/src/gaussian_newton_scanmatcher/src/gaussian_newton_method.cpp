#include <map.h>
#include "gaussian_newton_method.h"

const double GN_PI = 3.1415926;

//进行角度正则化．
double GN_NormalizationAngle(double angle)
{
    if(angle > GN_PI)
        angle -= 2*GN_PI;
    else if(angle < -GN_PI)
        angle += 2*GN_PI;

    return angle;
}

Eigen::Matrix3d GN_V2T(Eigen::Vector3d vec)
{
    Eigen::Matrix3d T;
    T  << cos(vec(2)),-sin(vec(2)),vec(0),
            sin(vec(2)), cos(vec(2)),vec(1),
            0,           0,     1;

    return T;
}

//对某一个点进行转换．
Eigen::Vector2d GN_TransPoint(Eigen::Vector2d pt,Eigen::Matrix3d T)
{
    Eigen::Vector3d tmp_pt(pt(0),pt(1),1);
    tmp_pt = T * tmp_pt;
    return Eigen::Vector2d(tmp_pt(0),tmp_pt(1));
}



//用激光雷达数据创建势场．
map_t* CreateMapFromLaserPoints(Eigen::Vector3d map_origin_pt,
                                std::vector<Eigen::Vector2d> laser_pts,
                                double resolution)
{
    map_t* map = map_alloc();

    map->origin_x = map_origin_pt(0);
    map->origin_y = map_origin_pt(1);
    map->resolution = resolution;

    //固定大小的地图，必要时可以扩大．
    map->size_x = 10000;
    map->size_y = 10000;

    map->cells = (map_cell_t*)malloc(sizeof(map_cell_t)*map->size_x*map->size_y);

    //高斯平滑的sigma－－固定死
    map->likelihood_sigma = 0.5;

    Eigen::Matrix3d Trans = GN_V2T(map_origin_pt);

    //设置障碍物
    for(int i = 0; i < laser_pts.size();i++)
    {
        Eigen::Vector2d tmp_pt = GN_TransPoint(laser_pts[i],Trans);

        int cell_x,cell_y;
        cell_x = MAP_GXWX(map,tmp_pt(0));
        cell_y = MAP_GYWY(map,tmp_pt(1));

        map->cells[MAP_INDEX(map,cell_x,cell_y)].occ_state = CELL_STATUS_OCC;
    }

    //进行障碍物的膨胀--最大距离固定死．
    map_update_cspace(map,0.5);

    return map;
}


/**
 * @brief InterpMapValueWithDerivatives
 * 在地图上的进行插值，得到coords处的势场值和对应的关于位置的梯度．
 * 返回值为Eigen::Vector3d ans
 * ans(0)表示市场值
 * ans(1:2)表示梯度
 * @param map
 * @param coords
 * @return
 */
Eigen::Vector3d InterpMapValueWithDerivatives(map_t* map,Eigen::Vector2d& coords)
{
    Eigen::Vector3d ans;
    //TODO
    int x0=map->size_x/2.0+floor((coords(0)-map->origin_x)/map->resolution);
    int y0=map->size_y/2.0+floor((coords(1)-map->origin_y)/map->resolution);

    double P00,P10,P01,P11;
    P00=map->cells[MAP_INDEX(map,x0,y0)].score;
    P10=map->cells[MAP_INDEX(map,x0+1,y0)].score;
    P01=map->cells[MAP_INDEX(map,x0,y0+1)].score;
    P11=map->cells[MAP_INDEX(map,x0+1,y0+1)].score;
    double u,v;
    u=(coords(0)-map->origin_x)/map->resolution+map->size_x/2.0-(double)x0;
    v=(coords(1)-map->origin_y)/map->resolution+map->size_y/2.0-(double)y0;

    ans(0)=P00*(1-u)*(1-v)+P10*u*(1-v)+P11*u*v+P01*(1-u)*v;
    ans(1)=(v*(P11-P01)+(1-v)*(P01-P00))/map->resolution;
    ans(2)=(u*(P11-P01)-(1-u)*(P01-P00))/map->resolution;
    //END OF TODO

    return ans;
}


/**
 * @brief ComputeCompleteHessianAndb
 * 计算H*dx = b中的H和b
 * @param map
 * @param now_pose
 * @param laser_pts
 * @param H
 * @param b
 */
void ComputeHessianAndb(map_t* map, Eigen::Vector3d now_pose,
                        std::vector<Eigen::Vector2d>& laser_pts,
                        Eigen::Matrix3d& H, Eigen::Vector3d& b)
{
    H = Eigen::Matrix3d::Zero();
    b = Eigen::Vector3d::Zero();
    //TODO
    int nLaserPtNum=laser_pts.size();
    for(int i=0;i<nLaserPtNum;i++)
    {
        //dM
        Eigen::Vector2d xyPt=laser_pts[i];
        Eigen::Matrix3d T=GN_V2T(now_pose);
        Eigen::Vector2d ST=GN_TransPoint(xyPt,T);
        Eigen::Vector3d MdM=InterpMapValueWithDerivatives(map,ST);
        //dS/dT
        Eigen::Matrix<double,2,3> dSdT;
        double t1,t2;
        t1=-sin(now_pose(2))*xyPt(0)-cos(now_pose(2))*xyPt(1);
        t2=cos(now_pose(2))*xyPt(0)-sin(now_pose(2))*xyPt(1);
        dSdT<<1.0,0.0,t1;0.0,1.0,t2;
        //sum
        Eigen::Matrix3d deltaH;
        Eigen::Vector3d deltab;
        Eigen::Vector2d dM;
        double Mst;
        dM<<MdM(1),MdM(2);
        Mst=MdM(0);
        deltaH=(dM.transpose()*dSdT).transpose()*(dM.transpose()*dSdT);
        H=H+deltaH;

        deltab=(dM.transpose()*dSdT).transpose()*(1-Mst);
        b=b+deltab;
    }
    //END OF TODO
}


/**
 * @brief GaussianNewtonOptimization
 * 进行高斯牛顿优化．
 * @param map
 * @param init_pose
 * @param laser_pts
 */
void GaussianNewtonOptimization(map_t*map,Eigen::Vector3d& init_pose,std::vector<Eigen::Vector2d>& laser_pts)
{
    int maxIteration = 10;
    Eigen::Vector3d now_pose = init_pose;

    for(int i = 0; i < maxIteration;i++)
    {
        //TODO
        //compute H and b
        Eigen::Matrix3d H;
        Eigen::Vector3d b;
        ComputeHessianAndb(map,now_pose,laser_pts,H,b);
        Eigen::Vector3d deltaT = Eigen::Vector3d::Zero();
        //deltaT = H.inverse()*b;
        deltaT = H.ldlt().solve(b)*0.01;
        now_pose = now_pose +deltaT;
        if(deltaT.norm()<0.0001)
        {
            break;
        }
        //END OF TODO
    }
    init_pose = now_pose;

}
