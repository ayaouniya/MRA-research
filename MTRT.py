#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 22:03:07 2022

@author: jiomer
"""
#%%
#import nis
import numpy as np
import matplotlib.pyplot as plt

#%%

def Gr(theta,phi):
    """
    If needed, we can usinga simulated antenna gain plot, now all directions is 1
    
    """
    return 1

def Bv(Tb,v):
    """
    Transform the radio object to intensity
    v in unit MHz
    """
    c = 2.998e8
    v = v*10**6
    k = 1.38e-23
    
    return 2*k*v**2*Tb/c**2

def QAQ(theta_q1,phi_q1,theta_f,phi_f,theta_m,phi_m):
    theta_q = np.arccos(np.sin(theta_q1)*np.cos(phi_f-phi_q1)*np.sin(theta_f)+np.cos(theta_q1)*np.cos(theta_f))
    phi_q = np.arctan((np.sin(theta_q1)*np.cos(phi_f-phi_q1)*np.cos(theta_f)-np.cos(theta_q1)*np.sin(theta_f))/(np.sin(theta_q1)*np.sin(phi_f-phi_q1)))-np.arctan((-np.sin(theta_m)*np.sin(phi_m-phi_f)*np.cos(theta_f)+np.cos(theta_m)*np.sin(theta_f))/(np.sin(theta_m)*np.sin(phi_m-phi_f)))
    #phi_q = np.abs(phi_q)
    return theta_q,phi_q

def cof(theta_q,phi_q,lmd,R,q=10**18):
    R = R*1000
    re = 2.818e-15
    Gr = 1
    A = 1 #lmd**2/(4*np.pi)
    k = Gr*A*(lmd/(2*R*np.sin(theta_q)**2*np.sqrt(np.e/(re*q))+np.sqrt(lmd/np.pi)*np.sin(theta_q)*np.cos(phi_q))) #k是反射系数
    return k


def MTvec2ang(x,y,z):
    """
    把流星尾迹起点终点坐标转换为xyz坐标
    """
    l = np.sqrt(x**2 + y**2 + z**2)
    x = x/l
    y = y/l
    z = z/l
    theta = np.arccos(z)
    sin_phi = y/np.sin(theta)
    cos_phi = x/np.sin(theta)
    if sin_phi >=0:
        phi = np.arccos(cos_phi)
    else:
        tem_phi = np.arccos(cos_phi)
        phi = 2*np.pi - tem_phi

    return theta,phi

def reflection(is_vector,n_vector):
    """
    计算反射向量
    """
    is_vector = np.array(is_vector)/np.dot(np.array(is_vector),np.array(is_vector))**0.5
    n_vector = np.array(n_vector)/np.dot(np.array(n_vector),np.array(n_vector))**0.5
    reflected_v = is_vector - 2*(np.dot(is_vector,n_vector))*n_vector
    reflected_v = reflected_v/np.dot(reflected_v,reflected_v)**0.5*np.dot(is_vector,is_vector)**0.5
    return reflected_v

def ang2vec(theta,phi,rad=True):
    """
    输入角度转化为单位向量，theta为天顶距，phi为方位角。
    """
    if rad == True:
        x = np.sin(theta)*np.cos(phi)
        y = np.sin(theta)*np.sin(phi)
        z = np.cos(theta)
    else:
        theta = theta/180*np.pi
        phi = phi/180*np.pi
        x = np.sin(theta)*np.cos(phi)
        y = np.sin(theta)*np.sin(phi)
        z = np.cos(theta)
    return x,y,z

def vec2ang(vec_x,vec_y,vec_z):
    """
    向量转化为角度theta，phi。
    """
    # normalized vector
    l = np.sqrt(vec_x**2+vec_y**2+vec_z**2)
    x = vec_x/l
    y = vec_y/l
    z = vec_z/l
    if z == 1:
        theta = 0
        phi = 0
    else:
        theta = np.arccos(z)
        sin_theta = np.sqrt(1-z**2)
        #print("sin theta",sin_theta)
        cos_phi = x/sin_theta
        if cos_phi > 1:
            cos_phi = 1
        #print("cos theta",cos_phi)
        sin_phi = y/sin_theta
        if sin_phi > 1:
            sin_phi = 1
        if sin_phi >=0:
            phi = np.arccos(cos_phi)
        else:
            tem_phi = np.arccos(cos_phi)
            phi = 2*np.pi - tem_phi
        
    return theta,phi
        
    
    return theta, 2*np.pi - phi

def numinrange(num,lower,upper):
    """
    判断一个数是否在区间[lower，upper]中。
    """
    if lower <= num <= upper:
        return True
    elif num < lower:
        return False
    elif num > upper:
        return False
    else:
        return False

def numnorange(num,lower,upper):
    """
    判断一个数是否不再区间[lower,upper]中
    """
    if lower <= num <= upper:
        return False
    elif num < lower:
        return True
    elif num > upper:
        return True
    else:
        return True


def getwl(Mm,wlo):
    """
    Mm 是流星z轴在ws坐标中的向量，
    wlo是流星尾迹原点在ws坐标中的表示。
    以流星z轴建立局部坐标系world local
    找到local x轴使得向量o-wlo在local坐标系的xoz平面。
    o是world space的原点【0,0,0】
    x1,y1,z1是流星尾迹终点
    wlx，wly，wlz是world local的基向量
    """
    wlo = np.array(wlo)
    Mm = np.array(Mm)

    x0 = wlo[0]
    y0 = wlo[1]
    z0 = wlo[2]
    x1 = x0 + Mm[0] 
    y1 = y0 + Mm[1]
    z1 = z0 + Mm[2]
    
    zp = (x0**2+y0**2+z0**2-x0*x1-y0*y1-z0*z1)*Mm[2]/np.dot(Mm,Mm) + z0
    xp = Mm[0]/Mm[2]*(zp-z0)+x0
    yp = Mm[1]/Mm[2]*(zp-z0)+y0

    wlx = np.array([-xp,-yp,-zp])/np.dot( np.array([-xp,-yp,-zp]),np.array([-xp,-yp,-zp]) )**0.5
    tem_wly = np.array([ Mm[2]*yp - Mm[1]*zp, Mm[0]*zp - Mm[2]*xp, Mm[1]*xp - Mm[0]*yp ])
    wly = tem_wly/np.dot(tem_wly,tem_wly)**0.5
    wlz = Mm/np.dot(Mm,Mm)**0.5

    return wlx,wly,wlz

def get_trans(is_vector,wlx,wly,wlz,ws2wl=True):
    """
    is_vector: 需要进行基变换的向量
    wlx，wly，wlz是local坐标的基
    ws2wl = True：生成变换矩阵并输出在local坐标系下的向量。
    ws2wl = False：生成变换矩阵并输出在space坐标系下的向量。
    """
    is_vector = np.array(is_vector)
    inv_tm = np.array([ [ wlx[0],wly[0],wlz[0] ],[ wlx[1],wly[1],wlz[1] ],[ wlx[2],wly[2],wlz[2] ] ])
    inv_tm = np.matrix(inv_tm)
    if ws2wl:
        tm = inv_tm.I
        ts_vec = np.dot(tm, is_vector)
        ts_vec = np.array(ts_vec)
        ts_v = np.array([ ts_vec[0][0],ts_vec[0][1],ts_vec[0][2] ])
    else:
        tm = inv_tm
        ts_vec = np.dot(tm, is_vector)
        ts_vec = np.array(ts_vec)
        ts_v = np.array([ ts_vec[0][0],ts_vec[0][1],ts_vec[0][2] ])

    return ts_v

#%%

"""
Class Rays

从观测点射出的光线，方向用theta，phi表示。其中theta是天顶距，范围理论上为 0～90度， 天顶为0度。phi为方位角，似乎大多数都选择正东方为0度，逆时针方向为正。

以观测点op为XYZ坐标系的原点，op指向天顶为Z轴，正东为Y轴，正南为X轴。

"""

class Ray():
    """
    从world space坐标系的原点（观测点）发射出的光线。

    """
    def __init__(self,theta,phi):
        
        # op is ibservation position
        # self.op = [0,0,0]
        
        self.theta = theta/180*np.pi

        self.phi = phi/180*np.pi
        
        self.vectorlength = 1
        
        self.vectors = np.array([ np.sin(self.theta)*np.cos(self.phi),
        np.sin(self.theta)*np.sin(self.phi), np.cos(self.theta) ])
    """
    def __str__(self):
        
        return f"The gain for this ray is {Ray.ray_gain(self)}"
    """
    
    def rayvectors(self):
        """
        得到光线的方向向量
        """
        return self.vectors.reshape(3)*Ray.ray_gain(self)
    
    
    def ray_gain(self):
        """
        return a ray with length equals to gain value
        也许不需要这个东西，做天空蒙板的时候再说
        """
        gain = 1
        return self.vectorlength*gain
    def plot_ray(self):
        """
        没什么用，就是一开始计算以及画图的时候坐标系没有选好搞得头疼用这个来画图检查位置对不对。

        """
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111)
        cir1 = plt.Circle(xy=(0.0,0.0),radius=1,alpha=0.2)
        plt.scatter(0,0,marker="*",color="r",s=200,label="zenith")
        plt.scatter(0,1,marker="2",color="r",s=200,label="North")
        plt.scatter(0,-1,marker="1",color="g",s=200,label="South")
        plt.scatter(1,0,marker="4",color="r",s=200,label="East")
        plt.scatter(-1,0,marker="3",color="g",s=200,label="West")

        plt.scatter([self.vectors[0]],[self.vectors[1]],s=50,label="Ray directions at sky")

        plt.xlim(-1.05,1.05)
        plt.ylim(-1.05,1.05)
        plt.xticks([])
        plt.yticks([])
        ax.add_patch(cir1)
        plt.legend(fontsize=15,loc='best')
        plt.title("Ray in the sky",fontsize=26)


#%%
class MeteorTrail():
    """
    This class for Meteor trail
    
    """
    
    def __init__(self,position,mt_theta,mt_phi,length=100,nu=50e6):
        """
        position：start point, lower side circle center 在space坐标
        mt_theta,mt_phi是流星z轴的方向向量
        流星尾迹长度默认100km，我也不知道对不对反正看起来很长就对了。
        流星尾迹半径由notes p17的公式计算得到，根频率有关。
        """
        self.origin = np.array(position)
        self.theta = mt_theta/180*np.pi
        self.phi = mt_phi/180*np.pi
        self.length = length
        self.radius = np.sqrt( 2.8179e-15*1e18*(2.998e8/nu)**2/(np.pi**2*np.e) )/1000
        # 流星尾迹的方向向量
        self.nvector = np.array([ np.sin(self.theta)*np.cos(self.phi),
        np.sin(self.theta)*np.sin(self.phi),np.cos(self.theta) ])
        # meteor end point, upper side circle center
        self.end = self.origin + self.length*self.nvector
        
        
    def __str__(self):
        return f"The meteor trail you defined lowwer crculer center located at {self.origin}, normalized direction vector is {self.nvector}."
    
    def trans_origin(self):
        """
        将ws坐标的原点变换到wold local坐标系
        """
        wl_origin = -MeteorTrail.trans_coords(self,self.origin)
        return wl_origin

    def get_zvector(self):
        """
        得到流星在ws坐标中的Mm向量（终点-起点）
        """
        return self.nvector*self.length

    def trans_coords(self,vector,wswl=True):
        if wswl == True:
            localx,localy,localz = getwl(self.nvector*self.length,self.origin)
            trans_vector = get_trans(vector,localx,localy,localz)
            return trans_vector
        else:
            localx,localy,localz = getwl(self.nvector*self.length,self.origin)
            trans_vector = get_trans(vector,localx,localy,localz,ws2wl=False)
            return trans_vector
    def inscylinder(self,isv,isvp):
        """
        判断相交发生在local坐标系中
        isv: insert vector under world local coords
        isvp: insert vector start point under world local coords
        解一元二次方程看是否有解，有解的话可能相交
        """
        vx = isv[0]
        vy = isv[1]
        vz = isv[2]
        px = isvp[0]
        py = isvp[1]
        pz = isvp[2]
        a = vx**2+vy**2
        b = 2*( px*vx+py*vy )
        c = px**2+py**2 - self.radius**2
        delta = b**2 - 4*a*c
        if delta < 0:
            return False,[],[]
        else:
            t1 = (-b - np.sqrt(delta))/(2*a)
            t2 = (-b + np.sqrt(delta))/(2*a)
            tmin = min(t1,t2)
            tmax = max(t1,t2)
            instx1 = px+tmin*vx
            insty1 = py+tmin*vy
            instz1 = pz+tmin*vz
            instx2 = px+tmax*vx
            insty2 = py+tmax*vy
            instz2 = pz+tmax*vz

            return True,np.array([instx1,insty1,instz1]),np.array([instx2,insty2,instz2])

    def return_nvector(self,isv,isvp):
        """
        如果inscylinder返回有解的话判断是否在流星尾迹长度的范围内
        pz1是距离观测点较近的交点
        pz2是距离观测点较远的交点
        1 如果pz1 = 0 & pz2在范围内：相较于下底面，相交面法向量为【0,0,-1】
        2 如果pz1在范围内，pz2在下底面，相较于圆柱侧面，法向量为【insert_point_x,insert_point_y,0】
        3 如果pz1不再范围内，pz2在范围内，相交于下底面
        4 如果pz1<0 & pz2 > length，相交于下底面
        5 pz1在范围内，pz2在范围内，相交于圆柱侧面
        也许有几个重复了，但是能用就行.jpg
        """
        inst,instp1,instp2 = MeteorTrail.inscylinder(self,isv,isvp)
        if inst == False:
            return False,[]
        else:
            pz1 = instp1[2]
            pz2 = instp2[2]
            if (pz1 == 0) & numinrange(pz2,0,self.length):
                return True,np.array([0,0,-1])
            elif (pz2 == self.length) & (numinrange(pz1,0,self.length)):
                return True,np.array([instp2[0],instp2[1],0])/np.dot(np.array([instp2[0],instp2[1],0]),np.array([instp2[0],instp2[1],0]))**0.5
            elif numinrange(pz2,0,self.length) & numnorange(pz1,0,self.length):
                return True, np.array([0,0,-1])
            elif (pz1 < 0) & (pz2 > self.length):
                return True,np.array([0,0,-1])
            elif numinrange(pz1,0,self.length) & numinrange(pz2,0,self.length):
                return True,np.array([instp1[0],instp1[1],0])/np.dot(np.array([instp1[0],instp1[1],0]),np.array([instp1[0],instp1[1],0]))**0.5
            else:
                return False,[]

# %%
