"""Analytic tests of the compiled slab formal solver (no kernel build)."""
import ctypes as ct
from pathlib import Path
import numpy as np
from polarized import build

def test():
    n,nx,na=51,140,5;kk=nx*na
    lib=ct.CDLL(str(build('build/transport', (n,nx,na))))
    fre=np.ctypeslib.as_array((ct.c_double*(4*nx+nx*(nx+1)//2)).in_dll(lib,'msqqwfre_'))
    opacity=np.geomspace(1e-4,1,nx)
    fre[3*nx+nx*(nx+1)//2:]=opacity
    cm=np.ctypeslib.as_array((ct.c_double*(nx+2*na+3*kk)).in_dll(lib,'msqqcm2_'))
    mu=(np.polynomial.legendre.leggauss(na)[0]+1)/2
    cm[nx:nx+na]=mu
    f=lib.msslabformal_;ptr=ct.POINTER(ct.c_double)
    f.argtypes=[ptr,ptr,ptr,ptr,ptr,ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    f.restype=None
    kappa=np.repeat(opacity,na)[:,None];cosine=np.tile(mu,nx)[:,None]
    errors=[]
    for h in [1e-6,.01,.1,1.]:
        s=np.ones((n,kk),order='F');out=np.zeros((kk,n),order='F');down=out.copy(order='F')
        hh=ct.c_double(h);order=ct.c_int(2);iq=ct.c_int(1)
        f(s.ctypes.data_as(ptr),s.ctypes.data_as(ptr),out.ctypes.data_as(ptr),down.ctypes.data_as(ptr),ct.byref(hh),ct.byref(order),ct.byref(iq))
        distance=h*np.arange(n)[None,:]
        exact=-np.expm1(-kappa*distance/cosine)/kappa
        err=np.max(abs(out-exact)/(1+exact));errors.append(float(err))
        assert err<2e-12,(h,err)
        assert np.allclose(down,out[:,::-1],rtol=2e-13,atol=2e-13)
        assert np.min(out)>=0
    # Unscattered Lambert injection must be exact even across an opaque cell.
    sw=np.ctypeslib.as_array((ct.c_int*7).in_dll(lib,'msqwitch_'));sw[6]=0
    seed=np.ctypeslib.as_array((ct.c_double*nx).in_dll(lib,'msqbb_'));seed[:]=1
    order=ct.c_int(1);iq=ct.c_int(0)
    f(s.ctypes.data_as(ptr),s.ctypes.data_as(ptr),out.ctypes.data_as(ptr),down.ctypes.data_as(ptr),ct.byref(hh),ct.byref(order),ct.byref(iq))
    assert np.allclose(out,np.exp(-kappa*distance/cosine),atol=2e-14,rtol=2e-13)
    assert np.max(abs(down))==0
    print('Compiled constant-source / symmetry / positivity tests passed:',errors)

def test_passivity():
    # Cold, scalar, isotropic limit isolates spatial quadrature from the
    # energy kernel and every polarization coupling. Vacuum boundaries.
    mu,w=np.polynomial.legendre.leggauss(9);mu=(mu+1)/2;w=w/2
    for tau in [.1,1.,3.,5.]:
        n=51;h=tau/(n-1);old=np.zeros((n,n));new=old.copy()
        for m,weight in zip(mu,w):
            d=h/m;r=np.exp(-d);wu=-np.expm1(-d)/d-r;wd=1+np.expm1(-d)/d
            up=np.zeros((n,n));down=up.copy()
            for j in range(1,n):
                z=np.arange(j+1);a=np.ones(j+1);a[[0,-1]]=.5
                old[j,z]+=weight*.5*h/m*a*np.exp(-(j-z)*d)
                z=np.arange(n-j-1,n);a=np.ones(j+1);a[[0,-1]]=.5
                old[n-j-1,z]+=weight*.5*h/m*a*np.exp(-(z-(n-j-1))*d)
                up[j]=r*up[j-1];up[j,j-1]+=wu;up[j,j]+=wd
                k=n-j-1;down[k]=r*down[k+1];down[k,k+1]+=wu;down[k,k]+=wd
            new+=weight*.5*(up+down)
        ro=float(np.max(np.linalg.eigvals(old).real));rn=float(np.max(np.linalg.eigvals(new).real))
        assert 0<rn<1
        if tau==5:assert ro>1
        print(f'tau={tau:g}: old spectral radius {ro:.10f}; corrected {rn:.10f}')

def test_source_product():
    lib=ct.CDLL(str(build('build/transport', (51,140,5))))
    rng=np.random.default_rng(917);kk=700;n=51
    matrices=[np.asfortranarray(rng.normal(size=(kk,kk))/kk) for _ in range(4)]
    inputs=[np.asfortranarray(rng.normal(size=(kk,n))) for _ in range(4)]
    p,m,pq,mq=matrices;ip,im,qp,qm=inputs
    expected=[p.T@ip+m.T@im+pq.T@qp+mq.T@qm,m.T@ip+p.T@im+mq.T@qp+pq.T@qm]
    output=[np.zeros((kk,n),order='F') for _ in range(2)]
    ptr=ct.POINTER(ct.c_double);f=lib.mspolproduct_;f.argtypes=[ptr]*10;f.restype=None
    f(*[a.ctypes.data_as(ptr) for a in matrices+inputs+output])
    error=max(float(np.max(abs(a-b))) for a,b in zip(expected,output))
    assert error<1e-12,error
    print('Signed BLAS source product agrees with independent matrix expression:',error)

if __name__=='__main__':test();test_passivity();test_source_product()

