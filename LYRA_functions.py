import numpy as np
from scipy.spatial import ConvexHull
from sklearn.cluster import MeanShift, estimate_bandwidth
import pandas as pd
from matplotlib.path import Path
# import base64

################################################################ generic fuctions:
def is_number(n):
    try:
        float(n) 
    except ValueError:
        return False
    return True

def quadratic_sum(x,y):
	return np.sqrt(np.power(x,2) + np.power(y,2))

#annulus area:
def aarea(rmin, rmax):
	return np.pi*(np.power(rmax,2) - np.power(rmin,2))

def acircle(r):
	return np.pi*np.power(r,2)

#from barplot to histogram:
def bartohist(x,y):
    y_len = int(np.sum(np.around(y,decimals=0)))
    new_y = np.array(np.around(y,decimals=0), dtype=int)
    new_x = np.zeros(y_len)
    a=b=j=0
    index=True
    while(index==True):
        if a < new_y[j]:
            new_x[b] = x[j]
            b+=1
            a+=1
        else:
            j+=1
            a=0
        if b == len(new_x):
            index = False
    return new_x

def gaussian(x,mu,sigma,norm):
	return norm*np.exp(-(x-mu)**2/(2*np.power(sigma,2)))

# convex hull
def convex_hull(ra,dec):
    coords_input = np.stack((ra,dec), axis=1)
    hull_pts_input = ConvexHull(coords_input)
    points = coords_input[hull_pts_input.vertices]
    points = np.append(points, [points[0]], axis=0) #close the polygon
    return points[:,0], points[:,1]

# polygons area
def PolyArea(X, Y):
	X_ = X - np.mean(X)	
	Y_ = Y - np.mean(Y)
	correction = X_[-1]*Y_[0] - Y_[-1]*X_[0]
	main_area = np.dot(X_[:-1], Y_[1:]) - np.dot(Y_[:-1], X_[1:])
	corr_area = 0.5*np.abs(main_area+correction)
	return corr_area

# polylines intersection
def _rect_inter_inner(x1,x2):
    n1=x1.shape[0]-1
    n2=x2.shape[0]-1
    X1=np.c_[x1[:-1],x1[1:]]
    X2=np.c_[x2[:-1],x2[1:]]    
    S1=np.tile(X1.min(axis=1),(n2,1)).T
    S2=np.tile(X2.max(axis=1),(n1,1))
    S3=np.tile(X1.max(axis=1),(n2,1)).T
    S4=np.tile(X2.min(axis=1),(n1,1))
    return S1,S2,S3,S4

def _rectangle_intersection_(x1,y1,x2,y2):
    S1,S2,S3,S4=_rect_inter_inner(x1,x2)
    S5,S6,S7,S8=_rect_inter_inner(y1,y2)

    C1=np.less_equal(S1,S2)
    C2=np.greater_equal(S3,S4)
    C3=np.less_equal(S5,S6)
    C4=np.greater_equal(S7,S8)

    ii,jj=np.nonzero(C1 & C2 & C3 & C4)
    return ii,jj

def intersection(x1,y1,x2,y2):
    """
	INTERSECTIONS Intersections of curves.
	Computes the (x,y) locations where two curves intersect.  The curves
	can be broken with NaNs or have vertical segments.
	usage:
	x,y=intersection(x1,y1,x2,y2)
    Example:
    a, b = 1, 2
    phi = np.linspace(3, 10, 100)
    x1 = a*phi - b*np.sin(phi)
    y1 = a - b*np.cos(phi)
    x2=phi    
    y2=np.sin(phi)+2
    x,y=intersection(x1,y1,x2,y2)
    plt.plot(x1,y1,c='r')
    plt.plot(x2,y2,c='g')
    plt.plot(x,y,'*k')
    plt.show()
    """
    # more points -> interpolate
    '''x1_new = []
    y1_new = []
    x2_new = []
    y2_new = []
    for i in range(len(x1_)-1):
        f1 = interpolate.interp1d([x1_[i], x1_[i+1]], [y1_[i], y1_[i+1]], kind='linear')
        x1_tmp = np.linspace(x1_[i], x1_[i+1], 10)
        y1_tmp = f1(x1_tmp)
        x1_new = np.append(x1_new, x1_tmp, axis=0)
        y1_new= np.append(y1_new, y1_tmp, axis=0)
    for i in range(len(x2_)-1):
        f2 = interpolate.interp1d([x2_[i], x2_[i+1]], [y2_[i], y2_[i+1]], kind='linear')
        x2_tmp = np.linspace(x2_[i], x2_[i+1], 10)
        y2_tmp = f2(x2_tmp)
        x2_new = np.append(x2_new, x2_tmp, axis=0)
        y2_new = np.append(y2_new, y2_tmp, axis=0)
    
    x1=x1_new
    y1=y1_new
    x2=x2_new
    y2=y2_new'''

    ii,jj=_rectangle_intersection_(x1,y1,x2,y2)
    n=len(ii)

    dxy1=np.diff(np.c_[x1,y1],axis=0)
    dxy2=np.diff(np.c_[x2,y2],axis=0)

    T=np.zeros((4,n))
    AA=np.zeros((4,4,n))
    AA[0:2,2,:]=-1
    AA[2:4,3,:]=-1
    AA[0::2,0,:]=dxy1[ii,:].T
    AA[1::2,1,:]=dxy2[jj,:].T

    BB=np.zeros((4,n))
    BB[0,:]=-x1[ii].ravel()
    BB[1,:]=-x2[jj].ravel()
    BB[2,:]=-y1[ii].ravel()
    BB[3,:]=-y2[jj].ravel()

    for i in range(n):
        try:
            T[:,i]=np.linalg.solve(AA[:,:,i],BB[:,i])
        except:
            T[:,i]=np.NaN

    in_range = (T[0,:] >=0) & (T[1,:] >=0) & (T[0,:] <=1) & (T[1,:] <=1)

    xy0=T[2:,in_range]
    xy0=xy0.T
    return xy0[:,0],xy0[:,1]

def find_pgons(ra, dec):
    X = np.stack((ra, dec), axis=1)
    if len(X) < 1000:
        ms = MeanShift(bandwidth=None, bin_seeding=True)
    elif (len(X) >= 1000) and (len(X) < 50000):
        bandwidth = estimate_bandwidth(X, quantile=0.02, n_samples=500, n_jobs=-1)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    else:
        bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=500, n_jobs=-1)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)

    ms.fit(X)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    print("number of estimated clusters : %d" % n_clusters_)

    area_cluster = 0
    pgons = []
    if n_clusters_ > 0:
        for k in range(n_clusters_):
            my_members = labels == k
            #cluster_center = cluster_centers[k]
            if len(X[my_members, 0]) > 3:
                points_ch_ra, points_ch_dec = convex_hull(X[my_members, 0], X[my_members, 1])
                area_tmp = PolyArea(points_ch_ra, points_ch_dec)
                area_cluster = area_cluster + area_tmp
                poly = Path(np.stack((points_ch_ra, points_ch_dec), axis=1))
                pgons = np.append(pgons, poly)
    else:            
        print("Warning: not enough cluster detected, go to convex hull estimation.")
    
    return pgons, area_cluster  # all the detected clusters and the total area

def print_output(path, message):
    outinfo = open(path,"a")
    outinfo.write("\n" + message + "\n")
    outinfo.close()

def chi_sqrd(N, x, mu, sigma):
    chi2 = []
    for i in range(N):
        chi2 = np.append(chi2, (x[i]-mu)**2/(sigma[i]**2))
        print(chi2[i])
    
    return np.sum(chi2)


################################################################ LR functions:

#f(r,c): angular separation prob. distribution, where
#sigma_o: 1sigma position error in onir srcs; sigma_x: 1sigma position error in X srcs; r=sqrt(ra_off**2+dec_off**2)
def FR(sigma_out, sigma_in, deltaRA, deltaDEC):
	sigma_ = quadratic_sum(sigma_in, sigma_out)
	r = quadratic_sum(deltaRA, deltaDEC)
	fr = (1/(2*np.pi*np.power(sigma_,2)))*np.exp(-(np.power(r,2))/(2*np.power(sigma_,2)))
	return fr

#LR: likelihood ratio
def LR(fr,ratio):
	return ratio*fr

#Re: realibiliy
def Re(summLR, LR_single, bayes, Q=1):
    if bayes == True:
        Q=1
    Re = LR_single/(summLR+(1-Q))
    return Re


################################################################ Bayes functions:

# everything in radians, steradians. These are multiplication factors
arcsec_to_rad = np.pi/(180*3600)
deg_to_rad = np.pi/180
deg2_to_sterad = (np.pi/180)**2
deg2_to_arcsec2 = 3600**2

def B_pos(sigma_out, sigma_in, deltaRA, deltaDEC):    # bayes positional factor
    # radians
    sigma2_input = np.power(sigma_in*arcsec_to_rad, 2) 
    sigma2_output = np.power(sigma_out*arcsec_to_rad, 2)
    psi = quadratic_sum(deltaRA*arcsec_to_rad, deltaDEC*arcsec_to_rad)
    bpos = (2/(sigma2_input+sigma2_output))*np.exp(-(np.power(psi,2))/(2*(sigma2_input+sigma2_output)))
    #print(sigma2_input,sigma2_output,psi)
    # arcsec
    # sigma2_input = np.power(sigma_in, 2)
    # sigma2_output = np.power(sigma_out, 2)
    # psi = quadratic_sum(deltaRA, deltaDEC)
    # bpos = (2/(sigma2_input+sigma2_output))*np.exp(-(np.power(psi,2))/(2*(sigma2_input+sigma2_output)))
    return bpos

def PH(dens1, dens2, Q=1): # prior probability
    dens1_st = dens1/deg2_to_sterad
    dens2_st = dens2/deg2_to_sterad
    #area = area*deg2_to_sterad
    dens1_st = dens1_st*4*np.pi
    dens2_st = dens2_st*4*np.pi
    ph = (dens2_st*Q)/(dens2_st*dens1_st)
    return ph

def PHD(PH, Bpos):    # posterior probability from Budavari and Szalay 2008
    phd = 1 / (1+((1-PH)/(Bpos*PH)))
    return phd

def Pany(PHD_sum):
    ph0d = 1
    phd_new = np.append(PHD_sum, ph0d)
    phd_norm = phd_new/np.sum(phd_new)
    pany = 1 - (phd_norm[-1]/np.sum(phd_norm))
    return pany

# check for duplicates
def check_duplicates_ids(filename):
    df = pd.read_csv(filename)
    out_ = df[df['C_flag'] == 1]['OutID'].duplicated().any()
    if out_ == True:
        duplicates = df[df['C_flag'] == 1][df[df['C_flag'] == 1]['OutID'].duplicated(keep=False)]
        print("Duplicates:")
        print(duplicates)
    in_ = df[df['C_flag'] == 1]['InID'].duplicated().any()
    if in_ == True:
        duplicates = df[df['C_flag'] == 1][df[df['C_flag'] == 1]['InID'].duplicated(keep=False)]
        print("Duplicates:")
        print(duplicates)
    return out_, in_
def check_duplicates_radec(filename):
    df = pd.read_csv(filename)
    out_ = df[df['C_flag'] == 1][['RAout','DECout']].duplicated().any()
    if out_ == True:
        duplicates = df[df['C_flag'] == 1][df[df['C_flag'] == 1][['RAout','DECout']].duplicated(keep=False)]
        print("Duplicates:")
        print(duplicates)
    in_ = df[df['C_flag'] == 1][['RAin','DECin']].duplicated().any()
    if in_ == True:
        duplicates = df[df['C_flag'] == 1][df[df['C_flag'] == 1][['RAin','DECin']].duplicated(keep=False)]
        print("Duplicates:")
        print(duplicates)
    return out_, in_

def find_plateau(x, y, threshold):
    x = np.array(x)
    y = np.array(y)
    y = y[x<0.8]
    x = x[x<0.8]
    max_value = np.max(y)
    max_index = np.argmax(y)

    if max_index > 0:  # Evaluate slope at maximum position
        slopes = []
        
        # Evaluate change of slope starting from the maximum value and going backwards
        for i in range(max_index, 1, -1):
            if x[i] != x[max_index]:
                slope2 = (y[i-1] - y[max_index]) / (x[i-1] - x[max_index])
                slopes.append(slope2)
        
        slope_std = np.std(slopes)
        threshold_value = threshold * slope_std
        
        for i in range(max_index, 1, -1):
            if x[i] != x[max_index]:
                slope2 = (y[i-1] - y[max_index]) / (x[i-1] - x[max_index])
                slope_change = slope2

                if slope_change > threshold_value:
                    max_value = y[i-1]
                    max_index = i-1
                    break
                
    # If no significant change of slope is found, return the maximum value and corresponding index
    return max_value, max_index



################################################################ AGNDB source names:
def LYRA_source_names(ra_deg, dec_deg):
    
    def deg2sex(ra='', dec='', round=False):
        RA, DEC, rs, ds = '', '', '', ''
        if dec:
            if str(dec)[0] == '-':
                ds, dec = '-', abs(dec)
            else:
                ds, dec = '+', abs(dec)
            deg = int(dec)
            decM = abs(int((dec-deg)*60))
            if round:
                decS = int((abs((dec-deg)*60)-decM)*60)
            else:
                decS = (abs((dec-deg)*60)-decM)*60
                decS = np.round(decS, 2)
            DEC = '{0}{1}{2}{3}'.format(ds, deg, decM, decS)

        if ra:
            if str(ra)[0] == '-':
                rs, ra = '-', abs(ra)
                print("\n### Warning: RA should be only positive.###\n")
                sys.exit()
            #else:
            #	rs, ra = '+', abs(ra)
            raH = int(ra/15)
            raM = int(((ra/15)-raH)*60)
            if round:
                raS = int(((((ra/15)-raH)*60)-raM)*60)
            else:
                raS = ((((ra/15)-raH)*60)-raM)*60
                raS = np.round(raS, 2)
            RA = '{0}{1}{2}{3}'.format(rs, raH, raM, raS)

        if ra and dec:
            return (RA, DEC)
        else:
            return RA or DEC

    names = np.empty(len(ra_deg), dtype=object)
    for i in range(len(ra_deg)):
        ra_conv, dec_conv = deg2sex(ra_deg[i], dec_deg[i])
        names[i] = "AGNDB J"+ra_conv+dec_conv
    return names
    # ra_conv, dec_conv = deg2sex(ra_deg, dec_deg)
    # name = "AGNDB J"+ra_conv+dec_conv
    # return name

################################################################



