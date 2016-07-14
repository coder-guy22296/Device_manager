import numpy as _np
from scipy.misc import factorial as _factorial
from scipy import optimize as _optimize
from scipy.ndimage import gaussian_filter as gf

nm = [(0,0), #1
      (1,1),
      (1,-1),
      (2,0),
      (2,-2),
      (2,2),
      (3,-1),
      (3,1),
      (3,-3),
      (3,3),
      (4,0), #11
      (4,2),
      (4,-2),
      (4,4),
      (4,-4),
      (5,1),
      (5,-1),
      (5,3),
      (5,-3),
      (5,5),
      (5,-5), #21
      (6,0),
      (6,-2),
      (6,2),
      (6,-4),
      (6,4),
      (6,-6),
      (6,6),
      (7,-1),
      (7,1),
      (7,-3), #31
      (7,3),
      (7,-5),
      (7,5),
      (7,-7),
      (7,7),
      (8,0)]


def radial((n,m),rho):
    """
    Computes the radial part of the Zernike polynomials.

    Parameters
    ----------
    (n,m) : tuple of ints
        Indices n and m
    rho: float
        Radial distance 

    Returns
    -------

    """
    
    if n<0 or m<0:
        print 'Error: n,m have to be >0.'
        return 0

    summation = _np.zeros_like(rho)

    #If n-m is odd, value is 0
    if _np.mod((n-m),2):
        print "For n-m odd Zernike polynomials are zero."
        return summation
   
    if n == 0:
        return _np.ones_like(rho)
    elif n == 1:
        return rho
    elif n == 2:
        if m == 0:
            return 2*rho**2-1
        elif m == 2:
            return rho**2
    elif n == 3:
        if m == 1:
            return 3*rho**2-2*rho
        elif m == 3:
            return rho**3
    elif n == 4:
        if m == 0:
            return 6*rho**4-6*rho**2+1
        elif m == 2:
            return 4*rho**4-3*rho**2
        elif m == 4:
            return rho**4
    
    for k in range(0, ((n-m)/2)+1):
        temp_top = ((-1)**k) * _factorial(n-k)
        temp_bot = _factorial(k) * _factorial((n+m)/2. - k) * _factorial((n-m)/2. - k)
        rho_term = pow(rho, n - (2*k))
        summation = summation + (rho_term * (temp_top/temp_bot))

        return summation


def zernike((n,m),rho,phi):
    """
    Computes the Zernike polynomials.

    Parameters
    ----------
    (n,m) : tuple of ints
        Zernike indices n and m
    rho: float
        Radial distance
    phi: float
        Azimuthal angle (in radians)

    Returns
    -------

    """
    if n < abs(m):
        print 'Error: n has to be >= abs(m)!'
        return 0

    rad_part = _np.sqrt(n+1)*radial((n,abs(m)),rho)
    
    if m == 0:
        return rad_part
    elif m>0:
        return _np.sqrt(2) * rad_part * _np.cos(m*phi)
    elif m<0:
        return _np.sqrt(2) * rad_part * _np.sin(m*phi)

def gradient_poly(index, r, theta):
    if index==2:
        sx = zernike(nm[0], r, theta)
        sy = _np.zeros_like(sx)
    if index==3:
        sy = zernike(nm[0], r, theta)
        sx = _np.zeros_like(sy)
    if index==4:
        sx = (1./_np.sqrt(2))*zernike(nm[1],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[2],r,theta)
    if index==5:
        sx = (1./_np.sqrt(2))*zernike(nm[2],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[1],r,theta)
    if index==6:
        sx = (1./_np.sqrt(2))*zernike(nm[1],r,theta)
        sy = (-1./_np.sqrt(2))*zernike(nm[2],r,theta)
    if index==7:
        sx = 0.5*zernike(nm[4],r,theta)
        sy = 0.5*(_np.sqrt(2)*zernike(nm[3],r,theta) - zernike(nm[5],r,theta))
    if index==8:
        sx = 0.5*(_np.sqrt(2)*zernike(nm[3],r,theta) + zernike(nm[5],r,theta))
        sy = 0.5*zernike(nm[4],r,theta)
    if index==9:
        sx = (1./_np.sqrt(2))*zernike(nm[4],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[5],r,theta)
    if index==10:
        sx = (1./_np.sqrt(2))*zernike(nm[5],r,theta)
        sy = (-1./_np.sqrt(2))*zernike(nm[4],r,theta)
    if index==11:
        sx = (1./_np.sqrt(2))*zernike(nm[7],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[6],r,theta)
    if index==12:
        sx = 0.5*(zernike(nm[7],r,theta)+zernike(nm[9],r,theta))
        sy = 0.5*(-1*zernike(nm[6],r,theta)+zernike(nm[8],r,theta))
    if index==13:
        sx = 0.5*(zernike(nm[6],r,theta)+zernike(nm[8],r,theta))
        sy = 0.5*(zernike(nm[7],r,theta)-zernike(nm[9],r,theta))
    if index==14:
        sx = (1./_np.sqrt(2))*zernike(nm[9],r,theta)
        sy = (-1./_np.sqrt(2))*zernike(nm[8],r,theta)
    if index==15:
        sx = (1./_np.sqrt(2))*zernike(nm[8],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[9],r,theta)
    if index==16:
        sx = 0.5*(_np.sqrt(2)*zernike(nm[10],r,theta) + zernike(nm[11],r,theta))
        sy = 0.5*zernike(nm[12],r,theta)
    if index==17:
        sx = 0.5*zernike(nm[12],r,theta)
        sy = 0.5*(_np.sqrt(2)*zernike(nm[10],r,theta) - zernike(nm[11],r,theta))
    if index==18:
        sx = 0.5*(zernike(nm[11],r,theta)+zernike(nm[13],r,theta))
        sy = 0.5*(-1*zernike(nm[12],r,theta)+zernike(nm[14],r,theta))
    if index==19:
        sx = 0.5*(zernike(nm[12],r,theta)+zernike(nm[14],r,theta))
        sy = 0.5*(zernike(nm[11],r,theta)-zernike(nm[13],r,theta))
    if index==20:
        sx = (1./_np.sqrt(2))*zernike(nm[13],r,theta)
        sy = (-1./_np.sqrt(2))*zernike(nm[14],r,theta)
    if index==21:
        sx = (1./_np.sqrt(2))*zernike(nm[14],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[13],r,theta)
    if index==22:
        sx = (1./_np.sqrt(2))*zernike(nm[15],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[16],r,theta)
    if index==23:
        sx = 0.5*(zernike(nm[16],r,theta)+zernike(nm[18],r,theta))
        sy = 0.5*(zernike(nm[15],r,theta)-zernike(nm[17],r,theta))
    if index==24:
        sx = 0.5*(zernike(nm[15],r,theta)+zernike(nm[17],r,theta))
        sy = 0.5*(-1*zernike(nm[16],r,theta)+zernike(nm[18],r,theta))
    if index==25:
        sx = 0.5*(zernike(nm[18],r,theta)+zernike(nm[20],r,theta))
        sy = 0.5*(zernike(nm[17],r,theta)-zernike(nm[19],r,theta))
    if index==26:
        sx = 0.5*(zernike(nm[17],r,theta)+zernike(nm[19],r,theta))
        sy = 0.5*(-1*zernike(nm[18],r,theta)+zernike(nm[20],r,theta))
    if index==27:
        sx = (1./_np.sqrt(2))*zernike(nm[20],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[19],r,theta)
    if index==28:
        sx = (1./_np.sqrt(2))*zernike(nm[19],r,theta)
        sy = (-1./_np.sqrt(2))*zernike(nm[20],r,theta)
    if index==29:
        sx = 0.5*zernike(nm[22],r,theta)
        sy = 0.5*(_np.sqrt(2)*zernike(nm[21],r,theta) - zernike(nm[23],r,theta))
    if index==30:
        sx = 0.5*(_np.sqrt(2)*zernike(nm[21],r,theta) + zernike(nm[23],r,theta))
        sy = 0.5*zernike(nm[22],r,theta)
    if index==31:
        sx = 0.5*(zernike(nm[22],r,theta)+zernike(nm[24],r,theta))
        sy = 0.5*(zernike(nm[23],r,theta)-zernike(nm[25],r,theta))
    if index==32:
        sx = 0.5*(zernike(nm[23],r,theta)+zernike(nm[25],r,theta))
        sy = 0.5*(-1*zernike(nm[22],r,theta)+zernike(nm[24],r,theta))
    if index==33:
        sx = 0.5*(zernike(nm[24],r,theta)+zernike(nm[26],r,theta))
        sy = 0.5*(zernike(nm[25],r,theta)-zernike(nm[27],r,theta))
    if index==34:
        sx = 0.5*(zernike(nm[25],r,theta)+zernike(nm[27],r,theta))
        sy = 0.5*(-1*zernike(nm[24],r,theta)+zernike(nm[26],r,theta))
    if index==35:
        sx = (1./_np.sqrt(2))*zernike(nm[26],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[27],r,theta)
    if index==36:
        sx = (1./_np.sqrt(2))*zernike(nm[27],r,theta)
        sy = (-1./_np.sqrt(2))*zernike(nm[26],r,theta)
    if index==37:
        sx = (1./_np.sqrt(2))*zernike(nm[29],r,theta)
        sy = (1./_np.sqrt(2))*zernike(nm[28],r,theta)
    return _np.array([sx,sy])

def basic_set_gradient(p,r,theta):
    result = _np.zeros((2,r.shape[0],r.shape[1]))
    for i in range(len(p)):
        result += p[i]*gradient_poly(i+2,r,theta)
    return result
    '''
    return p[0]*gradient_poly(2,r,theta) + \
           p[1]*gradient_poly(3,r,theta) + \
           p[2]*gradient_poly(4,r,theta) + \
           p[3]*gradient_poly(5,r,theta) + \
           p[4]*gradient_poly(6,r,theta) + \
           p[5]*gradient_poly(7,r,theta) + \
           p[6]*gradient_poly(8,r,theta) + \
           p[7]*gradient_poly(9,r,theta) + \
           p[8]*gradient_poly(10,r,theta) + \
           p[9]*gradient_poly(11,r,theta) + \
           p[10]*gradient_poly(12,r,theta) + \
           p[11]*gradient_poly(13,r,theta) + \
           p[12]*gradient_poly(14,r,theta) + \
           p[13]*gradient_poly(15,r,theta) + \
           p[14]*gradient_poly(16,r,theta) + \
           p[15]*gradient_poly(17,r,theta)
    '''

def return_grad(data, filter_grad=False, grad_filter_max=0.5, gfsize=0):
    #r_dx = _np.zeros_like(r)
    #r_dx[1:,:] = r[1:,:]-r[:-1,:]
    #r_dy = _np.zeros_like(r)
    #r_dy[:,1:] = r[:,1:]-r[:,:-1]
    data_dx = _np.zeros_like(data)
    data_dx[1:,:] = -1*(data[1:,:]-data[:-1,:])
    data_dy = _np.zeros_like(data)
    data_dy[:,1:] = (data[:,1:]-data[:,:-1])
    data_dxdy = -1*_np.array([data_dx,data_dy])
    if filter_grad:
        wx = _np.where(abs(data_dxdy[0])>grad_filter_max)
        wy = _np.where(abs(data_dxdy[1])>grad_filter_max)
        #print "Length wx: ", len(wx[0])
        '''
        for i in range(len(wx[0])):
            wx0 = wx[0][i]
            wx1 = wx[1][i]
            surrounding = data_dxdy[0][wx0-1:wx0+2,wx1-1:wx1+2].mean()
            data_dxdy[0][wx0,wx1] = surrounding
        for i in range(len(wy[0])):
            wy0 = wy[0][i]
            wy1 = wy[1][i]
            surrounding = data_dxdy[1][wy0-1:wy0+2,wy1-1:wy1+2].mean()
            data_dxdy[1][wy0,wy1] = surrounding
        '''
        data_dxdy[0][wx[0],wx[1]] = 0
        data_dxdy[1][wy[0],wy[1]] = 0
    if gfsize>0:
        data_dxdy[0] = gf(data_dxdy[0],gfsize)
        data_dxdy[1] = gf(data_dxdy[1],gfsize)
    return data_dxdy

def decompose_gradient(num_modes,data,r,theta,filter_grad = False,grad_filter_max=0.5,gfsize=0):
    basis = []
    for i in range(num_modes):
        basis.append(gradient_poly(i+2,r,theta))
    '''
    basis = [gradient_poly(2,r,theta),
             gradient_poly(3,r,theta),
             gradient_poly(4,r,theta),
             gradient_poly(5,r,theta),
             gradient_poly(6,r,theta),
             gradient_poly(7,r,theta),
             gradient_poly(8,r,theta),
             gradient_poly(9,r,theta),
             gradient_poly(10,r,theta),
             gradient_poly(11,r,theta),
             gradient_poly(12,r,theta),
             gradient_poly(13,r,theta),
             gradient_poly(14,r,theta),
             gradient_poly(15,r,theta),
             gradient_poly(16,r,theta),
             gradient_poly(17,r,theta)]
    '''
    
    data_dxdy = return_grad(data,filter_grad,grad_filter_max,gfsize=gfsize)
    '''
    if filter_grad:
        wx = _np.where(abs(data_dxdy[0])>grad_filter_max)
        wy = _np.where(abs(data_dxdy[1])>grad_filter_max)
        data_dxdy[0][wx[0],wx[1]] = 0
        data_dxdy[1][wy[0],wy[1]] = 0
    '''
    #radius = _np.sum((r<1)[r.shape[0]/2,:])/2
    radius = 1
    mask = _np.array([r<1,r<1])
    den = _np.array([_np.sum(z*z*mask) for z in basis])
    num = _np.array([_np.sum(data_dxdy*radius*z*mask) for z in basis])
    result1 = num/den
    cov_mat = _np.array([[_np.sum(zi*zj*mask) for zi in basis] for zj in basis])
    cov_mat_inv = _np.linalg.pinv(cov_mat)
    inner_products = _np.array([_np.sum(data_dxdy*radius*zi*mask) for zi in basis])
    result2 = _np.dot(cov_mat_inv,inner_products)
    return result1,result2

def fit_to_basic_gradient_set(image,p0,r,theta,weight=None,filter_grad = False,**kwargs):
    data_dxdy = return_grad(image)
    if filter_grad:
        wx = _np.where(data_dxdy[0]>0.5)
        wy = _np.where(data_dxdy[1]>0.5)
        data_dxdy[0][wx[0],wx[1]] = 0
        data_dxdy[1][wy[0],wy[1]] = 0
    mask = _np.array([r<1,r<1])
    #radius = _np.sum((r<1)[r.shape[0]/2,:])/2
    #radius = 1
    data = mask * data_dxdy * radius
    sumsquares = []
    count = [0]
    def fitfunction(p):
        Z = basic_set_gradient(p,r,theta)
        return Z*mask
    def errorfunction(p):
        #print 'Iteration', count[0]
        errormap = fitfunction(p)[mask]-data[mask]
        if weight != None:
            errormap *= weight
        error = _np.ravel(errormap)
        sumsquares.append(_np.sum(error**2))
        #print 'Sum of squares:',sumsquares[-1]
        count[0] += 1
        return error
    p,success = _optimize.leastsq(errorfunction,p0,**kwargs)
    return p,success,sumsquares

def gradient_to_zern(grad_coeffs, radius):
    new_coeffs = _np.zeros_like(grad_coeffs)
    print "Length of grad_coeffs: ", len(grad_coeffs)
    for i,(n,m) in enumerate(nm):
        #print "i index: ", i
        if len(grad_coeffs)>i:
            if i==0:
                new_coeffs[i]=0
            elif n==m:
                new_coeffs[i] = grad_coeffs[i-1] / _np.sqrt(2*n*(n+1))
            else:
                temp = grad_coeffs[i-1] / _np.sqrt(4*n*(n+1))
                new_n = n+2
                w1 = _np.where(_np.array(nm)[:,0] == new_n)
                w2 = _np.where(_np.array(nm)[:,1] == m)
                w = _np.intersect1d(w1[0],w2[0])
                if len(w)>1:
                    print "Multiple intersections..."
                if len(w)==0:
                    #print "No hits"
                    new_coeffs[i]=0
                else:
                    index_prime = w[0]-1
                    print "prime index: ", index_prime
                    if len(grad_coeffs)<=index_prime:
                        new_coeffs[i] = 0
                    else:
                        temp2 = grad_coeffs[index_prime] / _np.sqrt(4*(n+1)*(n+2))
                        new_coeffs[i] = temp - temp2
    return new_coeffs*radius
                    

def basic_set(p,r,theta):
    result = _np.zeros_like(r)
    for i in range(0,len(p)):
        result += p[i]*zernike(nm[i],r,theta)
    return result
    '''
    return  p[0]*zernike(nm[0],r,theta) + \
            p[1]*zernike(nm[1],r,theta) + \
            p[2]*zernike(nm[2],r,theta) + \
            p[3]*zernike(nm[3],r,theta) + \
            p[4]*zernike(nm[4],r,theta) + \
            p[5]*zernike(nm[5],r,theta) + \
            p[6]*zernike(nm[6],r,theta) + \
            p[7]*zernike(nm[7],r,theta) + \
            p[8]*zernike(nm[8],r,theta) + \
            p[9]*zernike(nm[9],r,theta) + \
            p[10]*zernike(nm[10],r,theta) + \
            p[11]*zernike(nm[11],r,theta) + \
            p[12]*zernike(nm[12],r,theta) + \
            p[13]*zernike(nm[13],r,theta) + \
            p[14]*zernike(nm[14],r,theta)
    '''

def decompose(num_modes, data,r,theta):
    basis = []
    for i in range(num_modes):
        basis.append(zernike(nm[i],r,theta))
    '''
    basis = [zernike(nm[0],r,theta),
             zernike(nm[1],r,theta),
             zernike(nm[2],r,theta),
             zernike(nm[3],r,theta),
             zernike(nm[4],r,theta),
             zernike(nm[5],r,theta),
             zernike(nm[6],r,theta),
             zernike(nm[7],r,theta),
             zernike(nm[8],r,theta),
             zernike(nm[9],r,theta),
             zernike(nm[10],r,theta),
             zernike(nm[11],r,theta),
             zernike(nm[12],r,theta),
             zernike(nm[13],r,theta),
             zernike(nm[14],r,theta)]
    '''
    cov_mat = _np.array([[_np.sum(zi*zj*(r<1)) for zi in basis] for zj in basis])
    cov_mat_inv = _np.linalg.pinv(cov_mat)
    inner_products = _np.array([_np.sum(data*zi*(r<1)) for zi in basis])
    return _np.dot(cov_mat_inv,inner_products)

def decompose2(data,r,theta):
    basis = [zernike((0,0),r,theta),
             zernike((1,-1),r,theta),
             zernike((1,1),r,theta),
             zernike((2,0),r,theta),
             zernike((2,-2),r,theta),
             zernike((2,2),r,theta),
             zernike((3,-1),r,theta),
             zernike((3,1),r,theta),
             zernike((3,-3),r,theta),
             zernike((3,3),r,theta),
             zernike((4,0),r,theta),
             zernike((4,-2),r,theta),
             zernike((4,2),r,theta),
             zernike((4,-4),r,theta),
             zernike((4,4),r,theta)]
    den = _np.array([_np.sum(z*z*(r<1)) for z in basis])
    num = _np.array([_np.sum(data*z*(r<1)) for z in basis])
    return num/den

def fit_to_basic_set(image,p0,r,theta,weight=None,**kwargs):
    data = _np.copy(image)
    data[r>1] = 0
    sumsquares = []
    count = [0]
    def fitfunction(p):
        Z = basic_set(p,r,theta)
        Z[r>1] = 0
        return Z
    def errorfunction(p):
        #print 'Iteration', count[0]
        errormap = fitfunction(p)-data
        if weight != None:
            errormap *= weight
        error = _np.ravel(errormap)
        sumsquares.append(_np.sum(error**2))
        #print 'Sum of squares:',sumsquares[-1]
        count[0] += 1
        return error
    p,success = _optimize.leastsq(errorfunction,p0,**kwargs)
    return p,success,sumsquares
