#include "gaussian_newton.h"
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

#include <iostream>

const double GN_PI = 3.1415926;
double GN_NormalizationAngle(double angle)
{
    if(angle > GN_PI)
        angle -= 2*GN_PI;
    else if(angle < -GN_PI)
        angle += 2*GN_PI;

    return angle;
}

//位姿-->转换矩阵
Eigen::Matrix3d PoseToTrans(Eigen::Vector3d x)
{
    Eigen::Matrix3d trans;
    trans << cos(x(2)),-sin(x(2)),x(0),
             sin(x(2)), cos(x(2)),x(1),
                     0,         0,    1;

    return trans;
}


//转换矩阵－－＞位姿
Eigen::Vector3d TransToPose(Eigen::Matrix3d trans)
{
    Eigen::Vector3d pose;
    pose(0) = trans(0,2);
    pose(1) = trans(1,2);
    pose(2) = atan2(trans(1,0),trans(0,0));

    return pose;
}

//计算整个pose-graph的误差
double ComputeError(std::vector<Eigen::Vector3d>& Vertexs,
                    std::vector<Edge>& Edges)
{
    double sumError = 0;
    for(int i = 0; i < Edges.size();i++)
    {
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        Eigen::Matrix3d Xi = PoseToTrans(xi);
        Eigen::Matrix3d Xj = PoseToTrans(xj);
        Eigen::Matrix3d Z  = PoseToTrans(z);

        Eigen::Matrix3d Ei = Z.inverse() *  Xi.inverse() * Xj;

        Eigen::Vector3d ei = TransToPose(Ei);


        sumError += ei.transpose() * infoMatrix * ei;
    }
    return sumError;
}


/**
 * @brief CalcJacobianAndError
 *         计算jacobian矩阵和error
 * @param xi    fromIdx
 * @param xj    toIdx
 * @param z     观测值:xj相对于xi的坐标
 * @param ei    计算的误差
 * @param Ai    相对于xi的Jacobian矩阵
 * @param Bi    相对于xj的Jacobian矩阵
 */
void CalcJacobianAndError(Eigen::Vector3d xi,Eigen::Vector3d xj,Eigen::Vector3d z,
                          Eigen::Vector3d& ei,Eigen::Matrix3d& Ai,Eigen::Matrix3d& Bi)
{
    //TODO--Start
    Eigen::Matrix3d Xi=PoseToTrans(xi);
    Eigen::Matrix3d Xj=PoseToTrans(xj);
    Eigen::Matrix3d Xij=Xi.inverse()*Xj;
    Eigen::Matrix3d Tij=PoseToTrans(z);

    Eigen::Matrix2d Rij=Xij.block(0,0,2,2);
    Eigen::Matrix2d Ri=Xi.block(0,0,2,2);
    Eigen::Matrix2d Rj=Xj.block(0,0,2,2);
    Eigen::Vector2d tij=Tij.block(0,2,2,1);
    Eigen::Vector2d ti=xi.block(0,0,2,1);
    Eigen::Vector2d tj=xj.block(0,0,2,1);

    Bi.setZero();
    Bi.block(0,0,2,2)=Rij.transpose()*Ri.transpose();
    Bi(2,2)=1.0;

    Eigen::Matrix2d dRiDt=Eigen::Matrix2d::Zero();
    //double theta = atan2(Ri(1,0),Ri(0,0));
    dRiDt<<-sin(xi(2)),cos(xi(2)),-cos(xi(2)),-sin(xi(2));
    Ai.setZero();
    Ai.block(0,0,2,2)=-Rij.transpose()*Ri.transpose();
    Ai(2,2)=-1.0;
    Ai.block(0,2,2,1)=Rij.transpose()*dRiDt*(tj-ti);

    Eigen::Vector3d xij=TransToPose(Xij);
    ei.block(0,0,2,1)=Rij.transpose()*(Ri.transpose()*(tj-ti)-tij);
    ei(2)=xj(2)-xi(2)-z(2);
    ei(2)=GN_NormalizationAngle(ei(2));
    //TODO--end
}

/**
 * @brief LinearizeAndSolve
 *        高斯牛顿方法的一次迭代．
 * @param Vertexs   图中的所有节点
 * @param Edges     图中的所有边
 * @return          位姿的增量
 */
Eigen::VectorXd  LinearizeAndSolve(std::vector<Eigen::Vector3d>& Vertexs,
                                   std::vector<Edge>& Edges)
{
    //申请内存
    Eigen::MatrixXd H(Vertexs.size() * 3,Vertexs.size() * 3);
    Eigen::VectorXd b(Vertexs.size() * 3);

    H.setZero();
    b.setZero();

    //固定第一帧
    Eigen::Matrix3d I;
    I.setIdentity();
    H.block(0,0,3,3) += I;

    //构造H矩阵　＆ b向量
    for(int i = 0; i < Edges.size();i++)
    {
        //提取信息
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        //计算误差和对应的Jacobian
        Eigen::Vector3d ei;
        Eigen::Matrix3d Ai;
        Eigen::Matrix3d Bi;
        CalcJacobianAndError(xi,xj,z,ei,Ai,Bi);

        //TODO--Start
        Eigen::MatrixXd dH(Vertexs.size()*3,Vertexs.size()*3);
        Eigen::VectorXd db(Vertexs.size() * 3);

        dH.setZero();
        db.setZero();

        dH.block(3*tmpEdge.xi,3*tmpEdge.xi,3,3)=Ai.transpose()*infoMatrix*Ai;
        dH.block(3*tmpEdge.xi,3*tmpEdge.xj,3,3)=Ai.transpose()*infoMatrix*Bi;
        dH.block(3*tmpEdge.xj,3*tmpEdge.xi,3,3)=Bi.transpose()*infoMatrix*Ai;
        dH.block(3*tmpEdge.xj,3*tmpEdge.xj,3,3)=Bi.transpose()*infoMatrix*Bi;

        db.block(3*tmpEdge.xi,0,3,1)=Ai.transpose()*infoMatrix*ei;
        db.block(3*tmpEdge.xj,0,3,1)=Bi.transpose()*infoMatrix*ei;

        H+=dH;
        b+=db;
        //TODO--End
    }

    //求解
    Eigen::VectorXd dx;

    //TODO--Start
    dx=-H.colPivHouseholderQr().solve(b);
    //TODO-End

    return dx;
}











