c                          GLMnet (5/17/08)
c
c
c                 Elastic net with squared-error loss
c
c dense predictor matrix:
c
c call elnet(ka,parm,no,ni,x,y,w,njd,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,
c            lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file
c
c
c sparse predictor matrix:
c
c call spelnet(ka,parm,no,ni,x,ix,jx,y,w,njd,jd,vp,ne,nx,nlam,flmin,ulam,th
c             isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse column format
c
c
c other inputs:
c
c   ka = algorithm flag
c      ka=1 => covariance updating algorithm
c      ka=2 => naive algorithm
c   parm = family member index (0 <= parm <= 1)
c        = 0.0 => ridge
c        = 1.0 => lasso
c   no = number of observations
c   ni = number of predictor variables
c   y(no) = response vector
c   w(no)= observation weights
c   njd = number of variables to delete
c   jd(njd) = do not use variables jd(1)...jd(njd)
c   vp(ni) = relative penalties for each predictor variable
c      vp(j) = 0 => jth variable unpenalized
c   ne = maximum number of variables allowed to enter largest model
c        (stopping criterion)
c   nx = maximum number of variables allowed to enter all models
c        along path (memory allocation, nx > ne).
c   nlam = (maximum) number of lamda values
c   flmin = user control of lamda values (>=0)
c      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
c      flmin >= 1.0 => use supplied lamda values (see below)
c   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
c   thr = convergence threshold for each lamda solution.
c      iterations stop when the maximum standardized coefficient
c      change from the previous iteration is less than thr
c      (suggested value, thr=1.0e-4)
c   isd = standarization flag:
c      isd = 0 => regression on original predictor variables
c      isd = 1 => regression on standardized predictor variables
c      Note: output solutions always reference original
c            variables locations and scales.
c
c output:
c
c   lmu = actual number of lamda values (solutions)
c   a0(lmu) = intercept values for each solution
c   ca(nx,lmu) = compressed coefficient values for each solution
c   ia(nx) = pointers to compressed coefficients
c   nin(lmu) = number of compressed coefficients for each solution
c   rsq(lmu) = R**2 values for each solution
c   alm(lmu) = lamda values corresponding to each solution
c   nlp = total passes over the data summed over all lamda values
c   jerr = error flag:
c      jerr  = 0 => no error
c      jerr != 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 10000 => maxval(vp) <= 0.0
c
c
c Note: x, y and w are overwritten by programs
c
c
c least-squares utility routines:
c
c uncompress coefficient vector for particular solution:
c
c call uncomp(ni,ca,ia,nin,a)
c
c input:
c
c    ni = total number of predictor variables
c    ca(nx) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni) =  uncompressed coefficient vector
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call modval(a0,ca,ia,nin,n,x,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(n) = model predictions
c
c
c
c
c          Symmetric binomial/multinomial logistic elastic net
c
c
c dense predictor matrix:
c
c call lognet (parm,no,ni,nc,x,y,njd,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,
c              maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file
c
c
c sparse predictor matrix:
c
c call splognet (parm,no,ni,nc,x,ix,jx,y,njd,jd,vp,ne,nx,nlam,flmin,
c             ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse column format
c
c
c other inputs:
c
c   parm, no, ni, njd, jd, vp, ne, nx, nlam, flmin, ulam, thr, isd, same as
c
c   nc = number of classes (distinct outcome values)
c        nc=1 => binomial two-class logistic regression
c            (all output references class 1)
c   y(no,max(2,nc)) = number of each class at each design point(overwritten
c   maxit = maximum number of iterations allowed for any lamda value
c           (suggested value, maxit = 100)
c   kopt = optimization flag
c      kopt = 0 => Newton-Raphson
c      kpot = 1 => modified Newton-Raphson (recommended)
c
c
c output:
c
c   lmu, ia, nin, alm, nlp, same as above
c
c   a0(nc,lmu) = intercept values for each class at each solution
c   ca(nx,nc,lmu) = compressed coefficient values for each class at
c                each solution
c   dev(lmu) = fraction of explained devience for each solution
c   jerr = error flag
c      jerr = 0  => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8000 + k => null probability < 1.0e-5 for class k
c         jerr = 9000 + k => null probability for class k
c                            > 1.0 - 1.0e-5
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output returned
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations. Solutions for
c            larger lamdas returned
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value. Solutions for
c            larger lamdas returned
c
c
c
c logistic utilitity routines:
c
c uncompress coefficient vector for particular solution:
c
c call luncomp(ni,nx,nc,ca,ia,nin,a)
c
c input:
c
c    ni, nx, nc = same as above
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c    a(ni,nc) =  uncompressed coefficient vectors
c                 referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor vectors:
c
c call lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans);
c
c input:
c
c    nt = number of observations
c    x(nt,ni) = full (uncompressed) predictor vectors
c    nc, nx = same as above
c    a0(nc) = intercepts
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c ans(nc,nt) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    nc, nx = same as above
c    a0(nc) = intercept
c    ca(nx,nc) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(nc,n) = model predictions
c
c
c
c
      subroutine elnet  (ka,parm,no,ni,x,y,w,njd,jd,vp,ne,nx,nlam,flmin
     *,ulam,thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam)
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)
      integer jd(njd),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: vq;
      if(maxval(vp) .gt. 0.0)goto 10021
      jerr=10000
      return
10021 continue
      allocate(vq(1:ni),stat=jerr)
      if(jerr.ne.0) return
      vq=max(0.0,vp)
      vq=vq*ni/sum(vq)
      if(ka .ne. 1)goto 10041
      call elnetu  (parm,no,ni,x,y,w,njd,jd,vq,ne,nx,nlam,flmin,ulam,thr,
     *isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051
10041 continue
      call elnetn (parm,no,ni,x,y,w,njd,jd,vq,ne,nx,nlam,flmin,ulam,thr,
     *isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue
10031 continue
      deallocate(vq)
      return
      end
      subroutine elnetu  (parm,no,ni,x,y,w,njd,jd,vp,ne,nx,nlam,flmin,ulam,
     *thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)
      integer jd(njd),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam
      integer, dimension (:), allocatable :: ju
      allocate(g(1:ni),stat=jerr)
      allocate(xm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xs(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ju(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(vlam(1:nlam),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      call chkvars(no,ni,x,ju)
      if(njd.gt.0) ju(jd(1:njd))=0
      if(maxval(ju) .gt. 0)goto 10071
      jerr=7777
      return
10071 continue
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)
      if(jerr.ne.0) return
      if(flmin.ge.1.0) vlam=ulam/ys
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,thr,xv,  lm
     *u,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return
10080 do 10081 k=1,lmu
      alm(k)=ys*alm(k)
      nk=nin(k)
10090 do 10091 l=1,nk
      ca(l,k)=ys*ca(l,k)/xs(ia(l))
10091 continue
10092 continue
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))
10081 continue
10082 continue
      deallocate(xm,xs,g,ju,xv,vlam)
      return
      end
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)
      integer ju(ni)
      real, dimension (:), allocatable :: v
      allocate(v(1:no),stat=jerr)
      if(jerr.ne.0) return
      w=w/sum(w)
      v=sqrt(w)
10100 do 10101 j=1,ni
      if(ju(j).eq.0)goto 10101
      xm(j)=dot_product(w,x(:,j))
      x(:,j)=v*(x(:,j)-xm(j))
      xv(j)=dot_product(x(:,j),x(:,j))
10101 continue
10102 continue
      if(isd .ne. 0)goto 10121
      xs=1.0
      goto 10131
10121 continue
10140 do 10141 j=1,ni
      if(ju(j).eq.0)goto 10141
      xs(j)=sqrt(xv(j))
      x(:,j)=x(:,j)/xs(j)
10141 continue
10142 continue
      xv=1.0
10131 continue
10111 continue
      ym=dot_product(w,y)
      y=v*(y-ym)
      ys=sqrt(dot_product(y,y))
      y=y/ys
      g=0.0
10150 do 10151 j=1,ni
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))
10151 continue
10152 continue
      deallocate(v)
      return
      end
      subroutine elnet1 (beta,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,ulam,thr,
     *xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99
     *9)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)
      real, dimension (:), allocatable :: a,da
      integer, dimension (:), allocatable :: mm
      real, dimension (:,:), allocatable :: c
      allocate(c(1:ni,1:nx),stat=jerr)
      allocate(a(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(da(1:ni),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      bta=max(beta,1.0e-3)
      omb=1.0-bta
      if(flmin .ge. 1.0)goto 10171
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
10171 continue
      rsq=0.0
      a=0.0
      mm=0
      nlp=0
      nin=nlp
      iz=0
      mnl=min(mnlam,nlam)
10180 do 10181 m=1,nlam
      if(flmin .lt. 1.0)goto 10201
      alm=ulam(m)
      goto 10191
10201 if(m .le. 2)goto 10211
      alm=alm*alf
      goto 10191
10211 if(m .ne. 1)goto 10221
      alm=big
      goto 10231
10221 continue
      alm=0.0
10240 do 10241 j=1,ni
      if(ju(j).eq.0)goto 10241
      if(vp(j).le.0.0)goto 10241
      alm=max(alm,abs(g(j))/vp(j))
10241 continue
10242 continue
      alm=alf*alm/bta
10231 continue
10191 continue
      dem=alm*omb
      ab=alm*bta
      rsq0=rsq
      jz=1
10250 continue
10251 continue
      if(iz*jz.ne.0) go to 10260
      nlp=nlp+1
      dlx=0.0
10270 do 10271 k=1,ni
      if(ju(k).eq.0)goto 10271
      ak=a(k)
      u=g(k)+ak*xv(k)
      v=abs(u)-vp(k)*ab
      a(k)=0.0
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)
      if(a(k).eq.ak)goto 10271
      if(mm(k) .ne. 0)goto 10291
      nin=nin+1
      if(nin.gt.nx)goto 10272
10300 do 10301 j=1,ni
      if(ju(j).eq.0)goto 10301
      if(mm(j) .eq. 0)goto 10321
      c(j,nin)=c(k,mm(j))
      goto 10301
10321 continue
      if(j .ne. k)goto 10341
      c(j,nin)=xv(j)
      goto 10301
10341 continue
      c(j,nin)=dot_product(x(:,j),x(:,k))
10301 continue
10302 continue
      mm(k)=nin
      ia(nin)=k
10291 continue
      del=a(k)-ak
      rsq=rsq+del*(2.0*g(k)-del*xv(k))
      dlx=max(abs(del)/sqrt(xv(k)),dlx)
10350 do 10351 j=1,ni
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del
10351 continue
10352 continue
10271 continue
10272 continue
      if(dlx.lt.thr)goto 10252
      if(nin.gt.nx)goto 10252
10260 continue
      iz=1
      da(1:nin)=a(ia(1:nin))
10360 continue
10361 continue
      nlp=nlp+1
      dlx=0.0
10370 do 10371 l=1,nin
      k=ia(l)
      ak=a(k)
      u=g(k)+ak*xv(k)
      v=abs(u)-vp(k)*ab
      a(k)=0.0
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)
      if(a(k).eq.ak)goto 10371
      del=a(k)-ak
      rsq=rsq+del*(2.0*g(k)-del*xv(k))
      dlx=max(abs(del)/sqrt(xv(k)),dlx)
10380 do 10381 j=1,nin
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del
10381 continue
10382 continue
10371 continue
10372 continue
      if(dlx.lt.thr)goto 10362
      goto 10361
10362 continue
      da(1:nin)=a(ia(1:nin))-da(1:nin)
10390 do 10391 j=1,ni
      if(mm(j).ne.0)goto 10391
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))
10391 continue
10392 continue
      jz=0
      goto 10251
10252 continue
      if(nin.gt.nx)goto 10182
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))
      kin(m)=nin
      rsqo(m)=rsq
      almo(m)=alm
      lmu=m
      if(m.lt.mnl)goto 10181
      if(flmin.ge.1.0)goto 10181
      me=0
10400 do 10401 j=1,nin
      if(ao(j,m).ne.0.0) me=me+1
10401 continue
10402 continue
      if(me.gt.ne)goto 10182
      if(rsq-rsq0.lt.sml*rsq)goto 10182
      if(rsq.gt.rsqmax)goto 10182
10181 continue
10182 continue
      deallocate(a,mm,c,da)
      return
      end
      subroutine elnetn (parm,no,ni,x,y,w,njd,jd,vp,ne,nx,nlam,flmin,ulam,
     *thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam)
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)
      integer jd(njd),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: xm,xs,xv,vlam
      integer, dimension (:), allocatable :: ju
      allocate(xm(1:ni),stat=jerr)
      allocate(xs(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ju(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(vlam(1:nlam),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      call chkvars(no,ni,x,ju)
      if(njd.gt.0) ju(jd(1:njd))=0
      if(maxval(ju) .gt. 0)goto 10421
      jerr=7777
      return
10421 continue
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)
      if(jerr.ne.0) return
      if(flmin.ge.1.0) vlam=ulam/ys
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,thr,xv,lmu
     *,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return
10430 do 10431 k=1,lmu
      alm(k)=ys*alm(k)
      nk=nin(k)
10440 do 10441 l=1,nk
      ca(l,k)=ys*ca(l,k)/xs(ia(l))
10441 continue
10442 continue
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))
10431 continue
10432 continue
      deallocate(xm,xs,ju,xv,vlam)
      return
      end
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)
      integer ju(ni)
      real, dimension (:), allocatable :: v
      allocate(v(1:no),stat=jerr)
      if(jerr.ne.0) return
      w=w/sum(w)
      v=sqrt(w)
10450 do 10451 j=1,ni
      if(ju(j).eq.0)goto 10451
      xm(j)=dot_product(w,x(:,j))
      x(:,j)=v*(x(:,j)-xm(j))
      xv(j)=dot_product(x(:,j),x(:,j))
10451 continue
10452 continue
      if(isd .ne. 0)goto 10471
      xs=1.0
      goto 10481
10471 continue
10490 do 10491 j=1,ni
      if(ju(j).eq.0)goto 10491
      xs(j)=sqrt(xv(j))
      x(:,j)=x(:,j)/xs(j)
10491 continue
10492 continue
      xv=1.0
10481 continue
10461 continue
      ym=dot_product(w,y)
      y=v*(y-ym)
      ys=sqrt(dot_product(y,y))
      y=y/ys
      deallocate(v)
      return
      end
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,x
     *v,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)
      real, dimension (:), allocatable :: a
      integer, dimension (:), allocatable :: mm
      allocate(a(1:ni),stat=jerr)
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      bta=max(beta,1.0e-3)
      omb=1.0-bta
      if(flmin .ge. 1.0)goto 10511
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
10511 continue
      rsq=0.0
      a=0.0
      mm=0
      nlp=0
      nin=nlp
      iz=0
      mnl=min(mnlam,nlam)
10520 do 10521 m=1,nlam
      if(flmin .lt. 1.0)goto 10541
      alm=ulam(m)
      goto 10531
10541 if(m .le. 2)goto 10551
      alm=alm*alf
      goto 10531
10551 if(m .ne. 1)goto 10561
      alm=big
      goto 10571
10561 continue
      alm=0.0
10580 do 10581 j=1,ni
      if(ju(j).eq.0)goto 10581
      if(vp(j).le.0.0)goto 10581
      alm=max(alm,abs(dot_product(y,x(:,j)))/vp(j))
10581 continue
10582 continue
      alm=alf*alm/bta
10571 continue
10531 continue
      dem=alm*omb
      ab=alm*bta
      rsq0=rsq
      jz=1
10590 continue
10591 continue
      if(iz*jz.ne.0) go to 10260
      nlp=nlp+1
      dlx=0.0
10600 do 10601 k=1,ni
      if(ju(k).eq.0)goto 10601
      gk=dot_product(y,x(:,k))
      ak=a(k)
      u=gk+ak*xv(k)
      v=abs(u)-vp(k)*ab
      a(k)=0.0
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)
      if(a(k).eq.ak)goto 10601
      if(mm(k) .ne. 0)goto 10621
      nin=nin+1
      if(nin.gt.nx)goto 10602
      mm(k)=nin
      ia(nin)=k
10621 continue
      del=a(k)-ak
      rsq=rsq+del*(2.0*gk-del*xv(k))
      y=y-del*x(:,k)
      dlx=max(abs(del)/sqrt(xv(k)),dlx)
10601 continue
10602 continue
      if(dlx.lt.thr)goto 10592
      if(nin.gt.nx)goto 10592
10260 continue
      iz=1
10630 continue
10631 continue
      nlp=nlp+1
      dlx=0.0
10640 do 10641 l=1,nin
      k=ia(l)
      gk=dot_product(y,x(:,k))
      ak=a(k)
      u=gk+ak*xv(k)
      v=abs(u)-vp(k)*ab
      a(k)=0.0
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)
      if(a(k).eq.ak)goto 10641
      del=a(k)-ak
      rsq=rsq+del*(2.0*gk-del*xv(k))
      y=y-del*x(:,k)
      dlx=max(abs(del)/sqrt(xv(k)),dlx)
10641 continue
10642 continue
      if(dlx.lt.thr)goto 10632
      goto 10631
10632 continue
      jz=0
      goto 10591
10592 continue
      if(nin.gt.nx)goto 10522
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))
      kin(m)=nin
      rsqo(m)=rsq
      almo(m)=alm
      lmu=m
      if(m.lt.mnl)goto 10521
      if(flmin.ge.1.0)goto 10521
      me=0
10650 do 10651 j=1,nin
      if(ao(j,m).ne.0.0) me=me+1
10651 continue
10652 continue
      if(me.gt.ne)goto 10522
      if(rsq-rsq0.lt.sml*rsq)goto 10522
      if(rsq.gt.rsqmax)goto 10522
10521 continue
10522 continue
      deallocate(a,mm)
      return
      end
      subroutine chkvars(no,ni,x,ju)
      real x(no,ni)
      integer ju(ni)
10660 do 10661 j=1,ni
      ju(j)=0
      t=x(1,j)
10670 do 10671 i=2,no
      if(x(i,j).eq.t)goto 10671
      ju(j)=1
      goto 10672
10671 continue
10672 continue
10661 continue
10662 continue
      return
      end
      subroutine uncomp(ni,ca,ia,nin,a)
      real ca(*),a(ni)
      integer ia(*)
      a=0.0
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)
      return
      end
      subroutine modval(a0,ca,ia,nin,n,x,f)
      real ca(nin),x(n,*),f(n)
      integer ia(nin)
      f=a0
      if(nin.le.0) return
10680 do 10681 i=1,n
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))
10681 continue
10682 continue
      return
      end
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,njd,jd,vp,ne,nx,
     *nlam,flmin,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)
      integer ix(*),jx(*),jd(njd),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: vq;
      if(maxval(vp) .gt. 0.0)goto 10701
      jerr=10000
      return
10701 continue
      allocate(vq(1:ni),stat=jerr)
      if(jerr.ne.0) return
      vq=max(0.0,vp)
      vq=vq*ni/sum(vq)
      if(ka .ne. 1)goto 10721
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,njd,jd,vq,ne,nx,nlam,flmin,
     *ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10731
10721 continue
      call spelnetn (parm,no,ni,x,ix,jx,y,w,njd,jd,vq,ne,nx,nlam,flmin,
     *ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10731 continue
10711 continue
      deallocate(vq)
      return
      end
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,njd,jd,vp,ne,nx,nlam,
     *flmin,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)
      integer ix(*),jx(*),jd(njd),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam
      integer, dimension (:), allocatable :: ju
      allocate(g(1:ni),stat=jerr)
      allocate(xm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xs(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ju(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(vlam(1:nlam),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      call spchkvars(no,ni,x,ix,ju)
      if(njd.gt.0) ju(jd(1:njd))=0
      if(maxval(ju) .gt. 0)goto 10751
      jerr=7777
      return
10751 continue
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,jerr)
      if(jerr.ne.0) return
      if(flmin.ge.1.0) vlam=ulam/ys
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t
     *hr,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return
10760 do 10761 k=1,lmu
      alm(k)=ys*alm(k)
      nk=nin(k)
10770 do 10771 l=1,nk
      ca(l,k)=ys*ca(l,k)/xs(ia(l))
10771 continue
10772 continue
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))
10761 continue
10762 continue
      deallocate(xm,xs,g,ju,xv,vlam)
      return
      end
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,j
     *err)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)
      integer ix(*),jx(*),ju(ni)
      w=w/sum(w)
10780 do 10781 j=1,ni
      if(ju(j).eq.0)goto 10781
      jb=ix(j)
      je=ix(j+1)-1
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2
10781 continue
10782 continue
      if(isd .ne. 0)goto 10801
      xs=1.0
      goto 10811
10801 continue
10820 do 10821 j=1,ni
      if(ju(j).ne.0) xs(j)=sqrt(xv(j))
10821 continue
10822 continue
      xv=1.0
10811 continue
10791 continue
      ym=dot_product(w,y)
      y=y-ym
      ys=sqrt(dot_product(w,y**2))
      y=y/ys
      g=0.0
10830 do 10831 j=1,ni
      if(ju(j).eq.0)goto 10831
      jb=ix(j)
      je=ix(j+1)-1
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)
10831 continue
10832 continue
      return
      end
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99
     *9)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)
      real, dimension (:), allocatable :: a,da
      integer, dimension (:), allocatable :: mm
      real, dimension (:,:), allocatable :: c
      allocate(c(1:ni,1:nx),stat=jerr)
      allocate(a(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(da(1:ni),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      bta=max(beta,1.0e-3)
      omb=1.0-bta
      if(flmin .ge. 1.0)goto 10851
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
10851 continue
      rsq=0.0
      a=0.0
      mm=0
      nlp=0
      nin=nlp
      iz=0
      mnl=min(mnlam,nlam)
10860 do 10861 m=1,nlam
      if(flmin .lt. 1.0)goto 10881
      alm=ulam(m)
      goto 10871
10881 if(m .le. 2)goto 10891
      alm=alm*alf
      goto 10871
10891 if(m .ne. 1)goto 10901
      alm=big
      goto 10911
10901 continue
      alm=0.0
10920 do 10921 j=1,ni
      if(ju(j).eq.0)goto 10921
      if(vp(j).le.0.0)goto 10921
      alm=max(alm,abs(g(j))/vp(j))
10921 continue
10922 continue
      alm=alf*alm/bta
10911 continue
10871 continue
      dem=alm*omb
      ab=alm*bta
      rsq0=rsq
      jz=1
10930 continue
10931 continue
      if(iz*jz.ne.0) go to 10260
      nlp=nlp+1
      dlx=0.0
10940 do 10941 k=1,ni
      if(ju(k).eq.0)goto 10941
      ak=a(k)
      u=g(k)+ak*xv(k)
      v=abs(u)-vp(k)*ab
      a(k)=0.0
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)
      if(a(k).eq.ak)goto 10941
      if(mm(k) .ne. 0)goto 10961
      nin=nin+1
      if(nin.gt.nx)goto 10942
10970 do 10971 j=1,ni
      if(ju(j).eq.0)goto 10971
      if(mm(j) .eq. 0)goto 10991
      c(j,nin)=c(k,mm(j))
      goto 10971
10991 continue
      if(j .ne. k)goto 11011
      c(j,nin)=xv(j)
      goto 10971
11011 continue
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))
10971 continue
10972 continue
      mm(k)=nin
      ia(nin)=k
10961 continue
      del=a(k)-ak
      rsq=rsq+del*(2.0*g(k)-del*xv(k))
      dlx=max(abs(del)/sqrt(xv(k)),dlx)
11020 do 11021 j=1,ni
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del
11021 continue
11022 continue
10941 continue
10942 continue
      if(dlx.lt.thr)goto 10932
      if(nin.gt.nx)goto 10932
10260 continue
      iz=1
      da(1:nin)=a(ia(1:nin))
11030 continue
11031 continue
      nlp=nlp+1
      dlx=0.0
11040 do 11041 l=1,nin
      k=ia(l)
      ak=a(k)
      u=g(k)+ak*xv(k)
      v=abs(u)-vp(k)*ab
      a(k)=0.0
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)
      if(a(k).eq.ak)goto 11041
      del=a(k)-ak
      rsq=rsq+del*(2.0*g(k)-del*xv(k))
      dlx=max(abs(del)/sqrt(xv(k)),dlx)
11050 do 11051 j=1,nin
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del
11051 continue
11052 continue
11041 continue
11042 continue
      if(dlx.lt.thr)goto 11032
      goto 11031
11032 continue
      da(1:nin)=a(ia(1:nin))-da(1:nin)
11060 do 11061 j=1,ni
      if(mm(j).ne.0)goto 11061
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))
11061 continue
11062 continue
      jz=0
      goto 10931
10932 continue
      if(nin.gt.nx)goto 10862
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))
      kin(m)=nin
      rsqo(m)=rsq
      almo(m)=alm
      lmu=m
      if(m.lt.mnl)goto 10861
      if(flmin.ge.1.0)goto 10861
      me=0
11070 do 11071 j=1,nin
      if(ao(j,m).ne.0.0) me=me+1
11071 continue
11072 continue
      if(me.gt.ne)goto 10862
      if(rsq-rsq0.lt.sml*rsq)goto 10862
      if(rsq.gt.rsqmax)goto 10862
10861 continue
10862 continue
      deallocate(a,mm,c,da)
      return
      end
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,njd,jd,vp,ne,nx,nlam,
     *flmin,ulam,  thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam)
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)
      integer ix(*),jx(*),jd(njd),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: xm,xs,xv,vlam
      integer, dimension (:), allocatable :: ju
      allocate(xm(1:ni),stat=jerr)
      allocate(xs(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ju(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(vlam(1:nlam),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      call spchkvars(no,ni,x,ix,ju)
      if(njd.gt.0) ju(jd(1:njd))=0
      if(maxval(ju) .gt. 0)goto 11091
      jerr=7777
      return
11091 continue
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,jerr)
      if(jerr.ne.0) return
      if(flmin.ge.1.0) vlam=ulam/ys
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t
     *hr,xm,xs,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return
11100 do 11101 k=1,lmu
      alm(k)=ys*alm(k)
      nk=nin(k)
11110 do 11111 l=1,nk
      ca(l,k)=ys*ca(l,k)/xs(ia(l))
11111 continue
11112 continue
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))
11101 continue
11102 continue
      deallocate(xm,xs,ju,xv,vlam)
      return
      end
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,je
     *rr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)
      integer ix(*),jx(*),ju(ni)
      w=w/sum(w)
11120 do 11121 j=1,ni
      if(ju(j).eq.0)goto 11121
      jb=ix(j)
      je=ix(j+1)-1
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2
11121 continue
11122 continue
      if(isd .ne. 0)goto 11141
      xs=1.0
      goto 11151
11141 continue
11160 do 11161 j=1,ni
      if(ju(j).ne.0) xs(j)=sqrt(xv(j))
11161 continue
11162 continue
      xv=1.0
11151 continue
11131 continue
      ym=dot_product(w,y)
      y=y-ym
      ys=sqrt(dot_product(w,y**2))
      y=y/ys
      return
      end
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99
     *9)
      real y(no),w(no),x(*),vp(ni),ulam(nlam)
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)
      real, dimension (:), allocatable :: a
      integer, dimension (:), allocatable :: mm
      allocate(a(1:ni),stat=jerr)
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      bta=max(beta,1.0e-3)
      omb=1.0-bta
      if(flmin .ge. 1.0)goto 11181
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
11181 continue
      rsq=0.0
      a=0.0
      mm=0
      o=0.0
      nlp=0
      nin=nlp
      iz=0
      mnl=min(mnlam,nlam)
11190 do 11191 m=1,nlam
      if(flmin .lt. 1.0)goto 11211
      alm=ulam(m)
      goto 11201
11211 if(m .le. 2)goto 11221
      alm=alm*alf
      goto 11201
11221 if(m .ne. 1)goto 11231
      alm=big
      goto 11241
11231 continue
      alm=0.0
11250 do 11251 j=1,ni
      if(ju(j).eq.0)goto 11251
      if(vp(j).le.0.0)goto 11251
      jb=ix(j)
      je=ix(j+1)-1
      alm=max(alm,abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))
     * /(vp(j)*xs(j))))
11251 continue
11252 continue
      alm=alf*alm/bta
11241 continue
11201 continue
      dem=alm*omb
      ab=alm*bta
      rsq0=rsq
      jz=1
11260 continue
11261 continue
      if(iz*jz.ne.0) go to 10260
      nlp=nlp+1
      dlx=0.0
11270 do 11271 k=1,ni
      if(ju(k).eq.0)goto 11271
      jb=ix(k)
      je=ix(k+1)-1
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)
      ak=a(k)
      u=gk+ak*xv(k)
      v=abs(u)-vp(k)*ab
      a(k)=0.0
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)
      if(a(k).eq.ak)goto 11271
      if(mm(k) .ne. 0)goto 11291
      nin=nin+1
      if(nin.gt.nx)goto 11272
      mm(k)=nin
      ia(nin)=k
11291 continue
      del=a(k)-ak
      rsq=rsq+del*(2.0*gk-del*xv(k))
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)
      o=o+del*xm(k)/xs(k)
      dlx=max(abs(del)/sqrt(xv(k)),dlx)
11271 continue
11272 continue
      if(dlx.lt.thr)goto 11262
      if(nin.gt.nx)goto 11262
10260 continue
      iz=1
11300 continue
11301 continue
      nlp=nlp+1
      dlx=0.0
11310 do 11311 l=1,nin
      k=ia(l)
      jb=ix(k)
      je=ix(k+1)-1
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)
      ak=a(k)
      u=gk+ak*xv(k)
      v=abs(u)-vp(k)*ab
      a(k)=0.0
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)
      if(a(k).eq.ak)goto 11311
      del=a(k)-ak
      rsq=rsq+del*(2.0*gk-del*xv(k))
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)
      o=o+del*xm(k)/xs(k)
      dlx=max(abs(del)/sqrt(xv(k)),dlx)
11311 continue
11312 continue
      if(dlx.lt.thr)goto 11302
      goto 11301
11302 continue
      jz=0
      goto 11261
11262 continue
      if(nin.gt.nx)goto 11192
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))
      kin(m)=nin
      rsqo(m)=rsq
      almo(m)=alm
      lmu=m
      if(m.lt.mnl)goto 11191
      if(flmin.ge.1.0)goto 11191
      me=0
11320 do 11321 j=1,nin
      if(ao(j,m).ne.0.0) me=me+1
11321 continue
11322 continue
      if(me.gt.ne)goto 11192
      if(rsq-rsq0.lt.sml*rsq)goto 11192
      if(rsq.gt.rsqmax)goto 11192
11191 continue
11192 continue
      deallocate(a,mm)
      return
      end
      subroutine spchkvars(no,ni,x,ix,ju)
      real x(*)
      integer ix(*),ju(ni)
11330 do 11331 j=1,ni
      ju(j)=0
      jb=ix(j)
      nj=ix(j+1)-jb
      if(nj.eq.0)goto 11331
      je=ix(j+1)-1
      if(nj .ge. no)goto 11351
11360 do 11361 i=jb,je
      if(x(i).eq.0.0)goto 11361
      ju(j)=1
      goto 11362
11361 continue
11362 continue
      goto 11371
11351 continue
      t=x(jb)
11380 do 11381 i=jb+1,je
      if(x(i).eq.t)goto 11381
      ju(j)=1
      goto 11382
11381 continue
11382 continue
11371 continue
11341 continue
11331 continue
11332 continue
      return
      end
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)
      real ca(*),x(*),f(n)
      integer ia(*),ix(*),jx(*)
      f=a0
11390 do 11391 j=1,nin
      k=ia(j)
      kb=ix(k)
      ke=ix(k+1)-1
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)
11391 continue
11392 continue
      return
      end
      function row_prod(i,j,ia,ja,ra,w)
      integer ia(*),ja(*)
      real ra(*),w(*)
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(
     *i),ia(j+1)-ia(j),w)
      return
      end
      function dot(x,y,mx,my,nx,ny,w)
      real x(*),y(*),w(*)
      integer mx(*),my(*)
      i=1
      j=i
      s=0.0
11400 continue
11401 continue
11410 continue
11411 if(mx(i).ge.my(j))goto 11412
      i=i+1
      if(i.gt.nx) go to 11420
      goto 11411
11412 continue
      if(mx(i).eq.my(j)) go to 11430
11440 continue
11441 if(my(j).ge.mx(i))goto 11442
      j=j+1
      if(j.gt.ny) go to 11420
      goto 11441
11442 continue
      if(mx(i).eq.my(j)) go to 11430
      goto 11401
11430 continue
      s=s+w(mx(i))*x(i)*y(j)
      i=i+1
      if(i.gt.nx)goto 11402
      j=j+1
      if(j.gt.ny)goto 11402
      goto 11401
11402 continue
11420 continue
      dot=s
      return
      end
      subroutine lognet (parm,no,ni,nc,x,y,njd,jd,vp,ne,nx,nlam,flmin,
     *ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),vp(ni),ulam(nlam)
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)
      integer jd(njd),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: xm,xs,ww,vq
      integer, dimension (:), allocatable :: ju
      if(maxval(vp) .gt. 0.0)goto 11461
      jerr=10000
      return
11461 continue
      allocate(ww(1:no),stat=jerr)
      allocate(ju(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(vq(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xm(1:ni),stat=ierr)
      jerr=jerr+ierr
      if(isd .le. 0)goto 11481
      allocate(xs(1:ni),stat=ierr)
      jerr=jerr+ierr
11481 continue
      if(jerr.ne.0) return
      call chkvars(no,ni,x,ju)
      if(njd.gt.0) ju(jd(1:njd))=0
      if(maxval(ju) .gt. 0)goto 11501
      jerr=7777
      return
11501 continue
      vq=max(0.0,vp)
      vq=vq*ni/sum(vq)
11510 do 11511 i=1,no
      ww(i)=sum(y(i,:))
      y(i,:)=y(i,:)/ww(i)
11511 continue
11512 continue
      ww=ww/sum(ww)
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)
      if(nc .ne. 1)goto 11531
      call lognet2n(parm,no,ni,x,y(:,1),ww,ju,vq,ne,nx,nlam,flmin,ulam,t
     *hr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
      goto 11541
11531 continue
      call lognetn(parm,no,ni,nc,x,y,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,
     *  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
11541 continue
11521 continue
      if(jerr.gt.0) return
11550 do 11551 k=1,lmu
      nk=nin(k)
11560 do 11561 ic=1,nc
      if(isd .le. 0)goto 11581
11590 do 11591 l=1,nk
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))
11591 continue
11592 continue
11581 continue
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))
11561 continue
11562 continue
11551 continue
11552 continue
      deallocate(ww,ju,vq,xm)
      if(isd.gt.0) deallocate(xs)
      return
      end
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)
      real x(no,ni),w(no),xm(ni),xs(ni)
      integer ju(ni)
11600 do 11601 j=1,ni
      if(ju(j).eq.0)goto 11601
      xm(j)=dot_product(w,x(:,j))
      x(1:no,j)=x(1:no,j)-xm(j)
      if(isd .le. 0)goto 11621
      xs(j)=sqrt(dot_product(w*x(:,j),x(:,j)))
      x(1:no,j)=x(1:no,j)/xs(j)
11621 continue
11601 continue
11602 continue
      return
      end
      subroutine lognet2n(parm,no,ni,x,y,w,ju,vp,ne,nx,nlam,flmin,ulam,s
     *hr,  isd,maxit,kopt,lmu,a0,a,m,kin,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=
     *5, devmax=0.999)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)
      integer ju(ni),m(nx),kin(nlam)
      real, dimension (:), allocatable :: b,bs,v,r,xv,q
      integer, dimension (:), allocatable :: mm
      allocate(b(0:ni),stat=jerr)
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(bs(0:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(r(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(v(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(q(1:no),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      fmax=log(1.0/pmin-1.0)
      fmin=-fmax
      vmin=(1.0+pmin)*pmin*(1.0-pmin)
      bta=max(parm,1.0e-3)
      omb=1.0-bta
      q0=dot_product(w,y)
      if(q0 .gt. pmin)goto 11641
      jerr=8001
      return
11641 continue
      if(q0 .lt. 1.0-pmin)goto 11661
      jerr=9001
      return
11661 continue
      vi=q0*(1.0-q0)
      b(1:ni)=0.0
      b(0)=log(q0/(1.0-q0))
      v=vi*w
      r=w*(y-q0)
      dev1=-(b(0)*q0+log(1.0-q0))
      q=q0
      if(isd .le. 0)goto 11681
      xv=0.25
      goto 11691
11681 continue
11700 do 11701 j=1,ni
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)
11701 continue
11702 continue
11691 continue
11671 continue
      xmz=vi
      bs=0.0
      if(flmin .ge. 1.0)goto 11721
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
11721 continue
      m=0
      mm=0
      nlp=0
      nin=nlp
      mnl=min(mnlam,nlam)
11730 do 11731 ilm=1,nlam
      if(flmin .lt. 1.0)goto 11751
      al=ulam(ilm)
      goto 11741
11751 if(ilm .le. 2)goto 11761
      al=al*alf
      goto 11741
11761 if(ilm .ne. 1)goto 11771
      al=big
      goto 11781
11771 continue
      al=0.0
11790 do 11791 j=1,ni
      if(ju(j).eq.0)goto 11791
      if(vp(j).le.0.0)goto 11791
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))
11791 continue
11792 continue
      al=alf*al/bta
11781 continue
11741 continue
      al2=al*omb
      al1=al*bta
      nit=0
11800 continue
11801 continue
      bs(0)=b(0)
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))
11810 continue
11811 continue
      nlp=nlp+1
      dlx=0.0
11820 do 11821 k=1,ni
      if(ju(k).eq.0)goto 11821
      bk=b(k)
      gk=dot_product(r,x(:,k))
      u=gk+xv(k)*b(k)
      au=abs(u)-vp(k)*al1
      if(au .gt. 0.0)goto 11841
      b(k)=0.0
      goto 11851
11841 continue
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)
11851 continue
11831 continue
      d=b(k)-bk
      if(abs(d).le.0.0)goto 11821
      dlx=max(dlx,abs(d))
      r=r-d*v*x(:,k)
      if(mm(k) .ne. 0)goto 11871
      nin=nin+1
      if(nin.gt.nx)goto 11822
      mm(k)=nin
      m(nin)=k
11871 continue
11821 continue
11822 continue
      if(nin.gt.nx)goto 11812
      d=sum(r)/xmz
      if(d .eq. 0.0)goto 11891
      b(0)=b(0)+d
      dlx=max(dlx,abs(d))
      r=r-d*v
11891 continue
      if(dlx.lt.shr)goto 11812
11900 continue
11901 continue
      nlp=nlp+1
      dlx=0.0
11910 do 11911 l=1,nin
      k=m(l)
      bk=b(k)
      gk=dot_product(r,x(:,k))
      u=gk+xv(k)*b(k)
      au=abs(u)-vp(k)*al1
      if(au .gt. 0.0)goto 11931
      b(k)=0.0
      goto 11941
11931 continue
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)
11941 continue
11921 continue
      d=b(k)-bk
      if(abs(d).le.0.0)goto 11911
      dlx=max(dlx,abs(d))
      r=r-d*v*x(:,k)
11911 continue
11912 continue
      d=sum(r)/xmz
      if(d .eq. 0.0)goto 11961
      b(0)=b(0)+d
      dlx=max(dlx,abs(d))
      r=r-d*v
11961 continue
      if(dlx.lt.shr)goto 11902
      goto 11901
11902 continue
      goto 11811
11812 continue
      if(nin.gt.nx)goto 11802
      if(abs(b(0)-bs(0)) .ge. shr)goto 11981
      ix=0
11990 do 11991 j=1,nin
      if(abs(b(m(j))-bs(m(j))).lt.shr)goto 11991
      ix=1
      goto 11992
11991 continue
11992 continue
      if(ix.eq.0)goto 11802
11981 continue
12000 do 12001 i=1,no
      fi=b(0)
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))
      if(fi .ge. fmin)goto 12021
      q(i)=0.0
      goto 12011
12021 if(fi .le. fmax)goto 12031
      q(i)=1.0
      goto 12041
12031 continue
      q(i)=1.0/(1.0+exp(-fi))
12041 continue
12011 continue
12001 continue
12002 continue
      v=w*q*(1.0-q)
      xmz=sum(v)
      if(xmz.le.vmin)goto 11802
      r=w*(y-q)
      if(kopt .ne. 0)goto 12061
12070 do 12071 j=1,nin
      xv(m(j))=dot_product(v,x(:,m(j))**2)
12071 continue
12072 continue
12061 continue
      nit=nit+1
      if(nit .le. maxit)goto 12091
      jerr=-ilm
      return
12091 continue
      goto 11801
11802 continue
      if(nin .le. nx)goto 12111
      jerr=-10000-ilm
      goto 11732
12111 continue
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))
      kin(ilm)=nin
      a0(ilm)=b(0)
      alm(ilm)=al
      lmu=ilm
      devi=dev2(no,w,y,q,pmin)
      dev(ilm)=(dev1-devi)/dev1
      if(xmz.le.vmin)goto 11732
      if(ilm.lt.mnl)goto 11731
      if(flmin.ge.1.0)goto 11731
      me=0
12120 do 12121 j=1,nin
      if(a(j,ilm).ne.0.0) me=me+1
12121 continue
12122 continue
      if(me.gt.ne)goto 11732
      if(dev(ilm).gt.devmax)goto 11732
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 11732
11731 continue
11732 continue
      deallocate(b,bs,v,r,xv,q,mm)
      return
      end
      function dev2(n,w,y,p,pmin)
      real w(n),y(n),p(n)
      pmax=1.0-pmin
      s=0.0
12130 do 12131 i=1,n
      pi=min(max(pmin,p(i)),pmax)
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))
12131 continue
12132 continue
      dev2=s
      return
      end
      subroutine lognetn(parm,no,ni,nc,x,y,w,ju,vp,ne,nx,nlam,flmin,ulam
     *,shr,  isd,maxit,kopt,lmu,a0,a,m,kin,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),w(no),vp(ni),ulam(nlam)
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)
      integer ju(ni),m(nx),kin(nlam)
      double precision, dimension (:,:), allocatable :: q
      double precision, dimension (:), allocatable :: sxp
      real, dimension (:), allocatable :: di,v,r
      real, dimension (:,:), allocatable :: b,bs,xv
      integer, dimension (:), allocatable :: mm,is
      allocate(b(0:ni,1:nc),stat=jerr)
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(r(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(v(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(is(1:max(nc,ni)),stat=ierr)
      jerr=jerr+ierr
      allocate(sxp(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(di(1:no),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      pmax=1.0-pmin
      emin=pmin/pmax
      emax=1.0/emin
      pfm=(1.0+pmin)*pmin
      pfx=(1.0-pmin)*pmax
      vmin=pfm*pmax
      bta=max(parm,1.0e-3)
      omb=1.0-bta
      dev1=0.0
12140 do 12141 ic=1,nc
      q0=dot_product(w,y(:,ic))
      if(q0 .gt. pmin)goto 12161
      jerr =8000+ic
      return
12161 continue
      if(q0 .lt. 1.0-pmin)goto 12181
      jerr =9000+ic
      return
12181 continue
      vi=q0*(1.0-q0)
      v=vi*w
      b(1:ni,ic)=0.0
      b(0,ic)=log(q0)
      dev1=dev1-q0*b(0,ic)
12141 continue
12142 continue
      if(isd .le. 0)goto 12201
      xv=0.25
      goto 12211
12201 continue
12220 do 12221 j=1,ni
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)
12221 continue
12222 continue
12211 continue
12191 continue
      b(0,:)=b(0,:)-sum(b(0,:))/nc
      sxp=0.0
12230 do 12231 ic=1,nc
      q(:,ic)=exp(b(0,ic))
      sxp=sxp+q(:,ic)
12231 continue
12232 continue
      if(flmin .ge. 1.0)goto 12251
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
12251 continue
      m=0
      mm=0
      nin=0
      nlp=0
      mnl=min(mnlam,nlam)
      bs=0.0
12260 do 12261 ilm=1,nlam
      if(flmin .lt. 1.0)goto 12281
      al=ulam(ilm)
      goto 12271
12281 if(ilm .le. 2)goto 12291
      al=al*alf
      goto 12271
12291 if(ilm .ne. 1)goto 12301
      al=big
      goto 12311
12301 continue
      al=0.0
12320 do 12321 ic=1,nc
      r=w*(y(:,ic)-q(:,ic)/sxp)
12330 do 12331 j=1,ni
      if(ju(j).eq.0)goto 12331
      if(vp(j).le.0.0)goto 12331
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))
12331 continue
12332 continue
12321 continue
12322 continue
      al=alf*al/bta
12311 continue
12271 continue
      al2=al*omb
      al1=al*bta
      nit=0
12340 continue
12341 continue
      ix=0
      jx=ix
      ig=0
12350 do 12351 ic=1,nc
      bs(0,ic)=b(0,ic)
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)
      xmz=0.0
12360 do 12361 i=1,no
      pic=q(i,ic)/sxp(i)
      if(pic .ge. pfm)goto 12381
      pic=0.0
      v(i)=0.0
      goto 12371
12381 if(pic .le. pfx)goto 12391
      pic=1.0
      v(i)=0.0
      goto 12401
12391 continue
      v(i)=w(i)*pic*(1.0-pic)
      xmz=xmz+v(i)
12401 continue
12371 continue
      r(i)=w(i)*(y(i,ic)-pic)
12361 continue
12362 continue
      if(xmz.le.vmin)goto 12351
      ig=1
      if(kopt .ne. 0)goto 12421
12430 do 12431 j=1,nin
      xv(m(j),ic)=dot_product(v,x(:,m(j))**2)
12431 continue
12432 continue
12421 continue
12440 continue
12441 continue
      nlp=nlp+1
      dlx=0.0
12450 do 12451 k=1,ni
      if(ju(k).eq.0)goto 12451
      bk=b(k,ic)
      gk=dot_product(r,x(:,k))
      u=gk+xv(k,ic)*b(k,ic)
      au=abs(u)-vp(k)*al1
      if(au .gt. 0.0)goto 12471
      b(k,ic)=0.0
      goto 12481
12471 continue
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)
12481 continue
12461 continue
      d=b(k,ic)-bk
      if(abs(d).le.0.0)goto 12451
      dlx=max(dlx,abs(d))
      r=r-d*v*x(:,k)
      if(mm(k) .ne. 0)goto 12501
      nin=nin+1
      if(nin .le. nx)goto 12521
      jx=1
      goto 12452
12521 continue
      mm(k)=nin
      m(nin)=k
12501 continue
12451 continue
12452 continue
      if(jx.gt.0)goto 12442
      d=sum(r)/xmz
      if(d .eq. 0.0)goto 12541
      b(0,ic)=b(0,ic)+d
      dlx=max(dlx,abs(d))
      r=r-d*v
12541 continue
      if(dlx.lt.shr)goto 12442
12550 continue
12551 continue
      nlp=nlp+1
      dlx=0.0
12560 do 12561 l=1,nin
      k=m(l)
      bk=b(k,ic)
      gk=dot_product(r,x(:,k))
      u=gk+xv(k,ic)*b(k,ic)
      au=abs(u)-vp(k)*al1
      if(au .gt. 0.0)goto 12581
      b(k,ic)=0.0
      goto 12591
12581 continue
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)
12591 continue
12571 continue
      d=b(k,ic)-bk
      if(abs(d).le.0.0)goto 12561
      dlx=max(dlx,abs(d))
      r=r-d*v*x(:,k)
12561 continue
12562 continue
      d=sum(r)/xmz
      if(d .eq. 0.0)goto 12611
      b(0,ic)=b(0,ic)+d
      dlx=max(dlx,abs(d))
      r=r-d*v
12611 continue
      if(dlx.lt.shr)goto 12552
      goto 12551
12552 continue
      goto 12441
12442 continue
      if(jx.gt.0)goto 12352
      if(abs(b(0,ic)-bs(0,ic)).gt.shr) ix=1
      if(ix .ne. 0)goto 12631
12640 do 12641 j=1,nin
      if(abs(b(m(j),ic)-bs(m(j),ic)) .le. shr)goto 12661
      ix=1
      goto 12642
12661 continue
12641 continue
12642 continue
12631 continue
12670 do 12671 i=1,no
      fi=b(0,ic)
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))
      fi=min(max(exmn,fi),exmx)
      sxp(i)=sxp(i)-q(i,ic)
      q(i,ic)=min(max(emin*sxp(i),exp(dble(fi))),emax*sxp(i))
      sxp(i)=sxp(i)+q(i,ic)
12671 continue
12672 continue
12351 continue
12352 continue
      if(jx.gt.0)goto 12342
      if(ix.eq.0)goto 12342
      if(ig.eq.0)goto 12342
      s=-sum(b(0,:))/nc
      b(0,:)=b(0,:)+s
      di=s
12680 do 12681 j=1,nin
      l=m(j)
      if(vp(l) .gt. 0.0)goto 12701
      s=sum(b(l,:))/nc
      goto 12711
12701 continue
      s=elc(parm,nc,b(l,:),is)
12711 continue
12691 continue
      b(l,:)=b(l,:)-s
      di=di-s*x(:,l)
12681 continue
12682 continue
      di=exp(di)
      sxp=sxp*di
12720 do 12721 ic=1,nc
      q(:,ic)=q(:,ic)*di
12721 continue
12722 continue
      nit=nit+1
      if(nit .le. maxit)goto 12741
      jerr=-ilm
      return
12741 continue
      goto 12341
12342 continue
      if(jx .le. 0)goto 12761
      jerr=-10000-ilm
      goto 12262
12761 continue
      devi=0.0
12770 do 12771 ic=1,nc
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)
      a0(ic,ilm)=b(0,ic)
12780 do 12781 i=1,no
      if(y(i,ic).le.0.0)goto 12781
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))
12781 continue
12782 continue
12771 continue
12772 continue
      kin(ilm)=nin
      alm(ilm)=al
      lmu=ilm
      dev(ilm)=(dev1-devi)/dev1
      if(ig.eq.0)goto 12262
      if(ilm.lt.mnl)goto 12261
      if(flmin.ge.1.0)goto 12261
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12262
      if(dev(ilm).gt.devmax)goto 12262
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12262
12261 continue
12262 continue
      deallocate(sxp,b,bs,v,r,xv,q,mm,is)
      return
      end
      function elc(parm,n,a,m)
      real a(n)
      integer m(n)
      fn=n
      am=sum(a)/fn
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 12801
      elc=am
      return
12801 continue
12810 do 12811 i=1,n
      m(i)=i
12811 continue
12812 continue
      call psort7(a,m,1,n)
      if(a(m(1)) .ne. a(m(n)))goto 12831
      elc=a(1)
      return
12831 continue
      if(mod(n,2) .ne. 1)goto 12851
      ad=a(m(n/2+1))
      goto 12861
12851 continue
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))
12861 continue
12841 continue
      if(parm .ne. 1.0)goto 12881
      elc=ad
      return
12881 continue
      b1=min(am,ad)
      b2=max(am,ad)
      k2=1
12890 continue
12891 if(a(m(k2)).gt.b1)goto 12892
      k2=k2+1
      goto 12891
12892 continue
      k1=k2-1
12900 continue
12901 if(a(m(k2)).ge.b2)goto 12902
      k2=k2+1
      goto 12901
12902 continue
      r=parm/((1.0-parm)*fn)
      is=0
      sm=n-2*(k1-1)
12910 do 12911 k=k1,k2-1
      sm=sm-2.0
      s=r*sm+am
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 12931
      is=k
      goto 12912
12931 continue
12911 continue
12912 continue
      if(is .eq. 0)goto 12951
      elc=s
      return
12951 continue
      r2=2.0*r
      s1=a(m(k1))
      am2=2.0*am
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)
      elc=s1
12960 do 12961 k=k1+1,k2
      s=a(m(k))
      if(s.eq.s1)goto 12961
      c=r2*sum(abs(a-s))+s*(s-am2)
      if(c .ge. cri)goto 12981
      cri=c
      elc=s
12981 continue
      s1=s
12961 continue
12962 continue
      return
      end
      function nintot(ni,nx,nc,a,m,nin,is)
      real a(nx,nc)
      integer m(nx),is(ni)
      is=0
      nintot=0
12990 do 12991 ic=1,nc
13000 do 13001 j=1,nin
      k=m(j)
      if(is(k).ne.0)goto 13001
      if(a(j,ic).eq.0.0)goto 13001
      is(k)=k
      nintot=nintot+1
13001 continue
13002 continue
12991 continue
12992 continue
      return
      end
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)
      real ca(nx,nc),a(ni,nc)
      integer ia(nx)
      a=0.0
13010 do 13011 ic=1,nc
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)
13011 continue
13012 continue
      return
      end
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)
      integer ia(nx)
13020 do 13021 i=1,nt
13030 do 13031 ic=1,nc
      ans(ic,i)=a0(ic)
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1
     *:nin)))
13031 continue
13032 continue
13021 continue
13022 continue
      return
      end
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,njd,jd,vp,ne,nx,nlam,
     *flmin,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
      real x(*),y(no,max(2,nc)),vp(ni),ulam(nlam)
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)
      integer ix(*),jx(*),jd(njd),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: xm,xs,ww,vq
      integer, dimension (:), allocatable :: ju
      if(maxval(vp) .gt. 0.0)goto 13051
      jerr=10000
      return
13051 continue
      allocate(ww(1:no),stat=jerr)
      allocate(ju(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(vq(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xs(1:ni),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      call spchkvars(no,ni,x,ix,ju)
      if(njd.gt.0) ju(jd(1:njd))=0
      if(maxval(ju) .gt. 0)goto 13071
      jerr=7777
      return
13071 continue
      vq=max(0.0,vp)
      vq=vq*ni/sum(vq)
13080 do 13081 i=1,no
      ww(i)=sum(y(i,:))
      y(i,:)=y(i,:)/ww(i)
13081 continue
13082 continue
      ww=ww/sum(ww)
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)
      if(nc .ne. 1)goto 13101
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),ww,ju,vq,ne,nx,nlam,flm
     *in,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev,alm,nlp,je
     *rr)
      goto 13111
13101 continue
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,ww,ju,vq,ne,nx,nlam,flmin,
     *ulam,thr,  isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
13111 continue
13091 continue
      if(jerr.gt.0) return
13120 do 13121 k=1,lmu
      nk=nin(k)
13130 do 13131 ic=1,nc
      if(isd .le. 0)goto 13151
13160 do 13161 l=1,nk
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))
13161 continue
13162 continue
13151 continue
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))
13131 continue
13132 continue
13121 continue
13122 continue
      deallocate(ww,ju,vq,xm,xs)
      return
      end
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)
      real x(*),w(no),xm(ni),xs(ni)
      integer ix(*),jx(*),ju(ni)
13170 do 13171 j=1,ni
      if(ju(j).eq.0)goto 13171
      jb=ix(j)
      je=ix(j+1)-1
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))
      if(isd.gt.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j
     *)**2)
13171 continue
13172 continue
      if(isd.eq.0) xs=1.0
      return
      end
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,w,ju,vp,ne,nx,nlam,
     *flmin,ulam,shr,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev,alm,nlp,jer
     *r)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=
     *5, devmax=0.999)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)
      real xb(ni),xs(ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q
      integer, dimension (:), allocatable :: mm
      allocate(b(0:ni),stat=jerr)
      allocate(xm(0:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(bs(0:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(q(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(r(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(v(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(sc(1:no),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      fmax=log(1.0/pmin-1.0)
      fmin=-fmax
      vmin=(1.0+pmin)*pmin*(1.0-pmin)
      bta=max(parm,1.0e-3)
      omb=1.0-bta
      q0=dot_product(w,y)
      if(q0 .gt. pmin)goto 13191
      jerr=8001
      return
13191 continue
      if(q0 .lt. 1.0-pmin)goto 13211
      jerr=9001
      return
13211 continue
      vi=q0*(1.0-q0)
      b(1:ni)=0.0
      b(0)=log(q0/(1.0-q0))
      v=vi*w
      r=w*(y-q0)
      dev1=-(b(0)*q0+log(1.0-q0))
      q=q0
      if(isd .le. 0)goto 13231
      xv=0.25
      goto 13241
13231 continue
13250 do 13251 j=1,ni
      if(ju(j).eq.0)goto 13251
      jb=ix(j)
      je=ix(j+1)-1
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)
13251 continue
13252 continue
13241 continue
13221 continue
      xm(0)=vi
      nlp=0
      nin=nlp
      if(flmin .ge. 1.0)goto 13271
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
13271 continue
      m=0
      mm=0
      nin=0
      o=0.0
      svr=o
      mnl=min(mnlam,nlam)
      bs=0.0
13280 do 13281 ilm=1,nlam
      if(flmin .lt. 1.0)goto 13301
      al=ulam(ilm)
      goto 13291
13301 if(ilm .le. 2)goto 13311
      al=al*alf
      goto 13291
13311 if(ilm .ne. 1)goto 13321
      al=big
      goto 13331
13321 continue
      al=0.0
13340 do 13341 j=1,ni
      if(ju(j).eq.0)goto 13341
      if(vp(j).le.0.0)goto 13341
      jb=ix(j)
      je=ix(j+1)-1
      jn=ix(j+1)-ix(j)
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o
      gj=dot_product(sc(1:jn),x(jb:je))
      gj=(gj-svr*xb(j))/xs(j)
      al=max(al,abs(gj)/vp(j))
13341 continue
13342 continue
      al=alf*al/bta
13331 continue
13291 continue
      al2=al*omb
      al1=al*bta
      nit=0
13350 continue
13351 continue
      bs(0)=b(0)
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))
13360 continue
13361 continue
      nlp=nlp+1
      dlx=0.0
13370 do 13371 k=1,ni
      if(ju(k).eq.0)goto 13371
      jb=ix(k)
      je=ix(k+1)-1
      jn=ix(k+1)-ix(k)
      bk=b(k)
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o
      gk=dot_product(sc(1:jn),x(jb:je))
      gk=(gk-svr*xb(k))/xs(k)
      u=gk+xv(k)*b(k)
      au=abs(u)-vp(k)*al1
      if(au .gt. 0.0)goto 13391
      b(k)=0.0
      goto 13401
13391 continue
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)
13401 continue
13381 continue
      d=b(k)-bk
      if(abs(d).le.0.0)goto 13371
      dlx=max(dlx,abs(d))
      if(mm(k) .ne. 0)goto 13421
      nin=nin+1
      if(nin.gt.nx)goto 13372
      mm(k)=nin
      m(nin)=k
      sc(1:jn)=v(jx(jb:je))
      xm(k)=dot_product(sc(1:jn),x(jb:je))
13421 continue
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)
      o=o+d*(xb(k)/xs(k))
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)
13371 continue
13372 continue
      if(nin.gt.nx)goto 13362
      d=svr/xm(0)
      if(d .eq. 0.0)goto 13441
      b(0)=b(0)+d
      dlx=max(dlx,abs(d))
      r=r-d*v
13441 continue
      svr=svr-d*xm(0)
      if(dlx.lt.shr)goto 13362
13450 continue
13451 continue
      nlp=nlp+1
      dlx=0.0
13460 do 13461 l=1,nin
      k=m(l)
      jb=ix(k)
      je=ix(k+1)-1
      jn=ix(k+1)-ix(k)
      bk=b(k)
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o
      gk=dot_product(sc(1:jn),x(jb:je))
      gk=(gk-svr*xb(k))/xs(k)
      u=gk+xv(k)*b(k)
      au=abs(u)-vp(k)*al1
      if(au .gt. 0.0)goto 13481
      b(k)=0.0
      goto 13491
13481 continue
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)
13491 continue
13471 continue
      d=b(k)-bk
      if(abs(d).le.0.0)goto 13461
      dlx=max(dlx,abs(d))
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)
      o=o+d*(xb(k)/xs(k))
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)
13461 continue
13462 continue
      d=svr/xm(0)
      if(d .eq. 0.0)goto 13511
      b(0)=b(0)+d
      dlx=max(dlx,abs(d))
      r=r-d*v
13511 continue
      svr=svr-d*xm(0)
      if(dlx.lt.shr)goto 13452
      goto 13451
13452 continue
      goto 13361
13362 continue
      if(nin.gt.nx)goto 13352
      if(abs(b(0)-bs(0)) .ge. shr)goto 13531
      kx=0
13540 do 13541 j=1,nin
      if(abs(b(m(j))-bs(m(j))).lt.shr)goto 13541
      kx=1
      goto 13542
13541 continue
13542 continue
      if(kx.eq.0)goto 13352
13531 continue
      sc=b(0)
      b0=0.0
13550 do 13551 j=1,nin
      l=m(j)
      jb=ix(l)
      je=ix(l+1)-1
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)
      b0=b0-b(l)*xb(l)/xs(l)
13551 continue
13552 continue
      sc=sc+b0
13560 do 13561 i=1,no
      fi=sc(i)
      if(fi .ge. fmin)goto 13581
      q(i)=0.0
      goto 13571
13581 if(fi .le. fmax)goto 13591
      q(i)=1.0
      goto 13601
13591 continue
      q(i)=1.0/(1.0+exp(-fi))
13601 continue
13571 continue
13561 continue
13562 continue
      v=w*q*(1.0-q)
      r=w*(y-q)
      xm(0)=sum(v)
      if(xm(0).lt.vmin)goto 13352
      svr=sum(r)
      o=0.0
13610 do 13611 l=1,nin
      j=m(l)
      jb=ix(j)
      je=ix(j+1)-1
      jn=ix(j+1)-ix(j)
      sc(1:jn)=v(jx(jb:je))
      xm(j)=dot_product(sc(1:jn),x(jb:je))
      if(kopt .ne. 0)goto 13631
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2
13631 continue
13611 continue
13612 continue
      nit=nit+1
      if(nit .le. maxit)goto 13651
      jerr=-ilm
      return
13651 continue
      goto 13351
13352 continue
      if(nin .le. nx)goto 13671
      jerr=-10000-ilm
      goto 13282
13671 continue
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))
      kin(ilm)=nin
      a0(ilm)=b(0)
      alm(ilm)=al
      lmu=ilm
      devi=dev2(no,w,y,q,pmin)
      dev(ilm)=(dev1-devi)/dev1
      if(ilm.lt.mnl)goto 13281
      if(flmin.ge.1.0)goto 13281
      me=0
13680 do 13681 j=1,nin
      if(a(j,ilm).ne.0.0) me=me+1
13681 continue
13682 continue
      if(me.gt.ne)goto 13282
      if(dev(ilm).gt.devmax)goto 13282
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13282
      if(xm(0).lt.vmin)goto 13282
13281 continue
13282 continue
      deallocate(xm,b,bs,v,r,sc,xv,q,mm)
      return
      end
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,w,ju,vp,ne,nx,nlam,f
     *lmin,ulam,  shr,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev,alm,nlp,je
     *rr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)
      double precision, dimension (:,:), allocatable :: q
      double precision, dimension (:), allocatable :: sxp
      real, dimension (:), allocatable :: sc,xm,v,r
      real, dimension (:,:), allocatable :: b,bs,xv
      integer, dimension (:), allocatable :: mm,is
      allocate(b(0:ni,1:nc),stat=jerr)
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(xm(0:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(r(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(v(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(is(1:max(nc,ni)),stat=ierr)
      jerr=jerr+ierr
      allocate(sxp(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(sc(1:no),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      pmax=1.0-pmin
      emin=pmin/pmax
      emax=1.0/emin
      pfm=(1.0+pmin)*pmin
      pfx=(1.0-pmin)*pmax
      vmin=pfm*pmax
      bta=max(parm,1.0e-3)
      omb=1.0-bta
      dev1=0.0
13690 do 13691 ic=1,nc
      q0=dot_product(w,y(:,ic))
      if(q0 .gt. pmin)goto 13711
      jerr =8000+ic
      return
13711 continue
      if(q0 .lt. 1.0-pmin)goto 13731
      jerr =9000+ic
      return
13731 continue
      vi=q0*(1.0-q0)
      v=vi*w
      b(1:ni,ic)=0.0
      b(0,ic)=log(q0)
      xm(0)=vi
      dev1=dev1-q0*b(0,ic)
13691 continue
13692 continue
      if(isd .le. 0)goto 13751
      xv=0.25
      goto 13761
13751 continue
13770 do 13771 j=1,ni
      if(ju(j).eq.0)goto 13771
      jb=ix(j)
      je=ix(j+1)-1
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)
13771 continue
13772 continue
13761 continue
13741 continue
      b(0,:)=b(0,:)-sum(b(0,:))/nc
      sxp=0.0
13780 do 13781 ic=1,nc
      q(:,ic)=exp(b(0,ic))
      sxp=sxp+q(:,ic)
13781 continue
13782 continue
      if(flmin .ge. 1.0)goto 13801
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
13801 continue
      m=0
      mm=0
      nin=0
      nlp=0
      mnl=min(mnlam,nlam)
      bs=0.0
      svr=0.0
      o=0.0
13810 do 13811 ilm=1,nlam
      if(flmin .lt. 1.0)goto 13831
      al=ulam(ilm)
      goto 13821
13831 if(ilm .le. 2)goto 13841
      al=al*alf
      goto 13821
13841 if(ilm .ne. 1)goto 13851
      al=big
      goto 13861
13851 continue
      al=0.0
13870 do 13871 ic=1,nc
      v=q(:,ic)/sxp
      r=w*(y(:,ic)-v)
      v=w*v*(1.0-v)
13880 do 13881 j=1,ni
      if(ju(j).eq.0)goto 13881
      if(vp(j).le.0.0)goto 13881
      jb=ix(j)
      je=ix(j+1)-1
      jn=ix(j+1)-ix(j)
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))
      gj=dot_product(sc(1:jn),x(jb:je))
      gj=(gj-svr*xb(j))/xs(j)
      al=max(al,abs(gj)/vp(j))
13881 continue
13882 continue
13871 continue
13872 continue
      al=alf*al/bta
13861 continue
13821 continue
      al2=al*omb
      al1=al*bta
      nit=0
13890 continue
13891 continue
      ixx=0
      jxx=ixx
      ig=0
13900 do 13901 ic=1,nc
      bs(0,ic)=b(0,ic)
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)
      xm(0)=0.0
      svr=0.0
      o=0.0
13910 do 13911 i=1,no
      pic=q(i,ic)/sxp(i)
      if(pic .ge. pfm)goto 13931
      pic=0.0
      v(i)=0.0
      goto 13921
13931 if(pic .le. pfx)goto 13941
      pic=1.0
      v(i)=0.0
      goto 13951
13941 continue
      v(i)=w(i)*pic*(1.0-pic)
      xm(0)=xm(0)+v(i)
13951 continue
13921 continue
      r(i)=w(i)*(y(i,ic)-pic)
      svr=svr+r(i)
13911 continue
13912 continue
      if(xm(0).le.vmin)goto 13901
      ig=1
13960 do 13961 l=1,nin
      j=m(l)
      jb=ix(j)
      je=ix(j+1)-1
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))
      if(kopt .ne. 0)goto 13981
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2
13981 continue
13961 continue
13962 continue
13990 continue
13991 continue
      nlp=nlp+1
      dlx=0.0
14000 do 14001 k=1,ni
      if(ju(k).eq.0)goto 14001
      jb=ix(k)
      je=ix(k+1)-1
      jn=ix(k+1)-ix(k)
      bk=b(k,ic)
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))
      gk=dot_product(sc(1:jn),x(jb:je))
      gk=(gk-svr*xb(k))/xs(k)
      u=gk+xv(k,ic)*b(k,ic)
      au=abs(u)-vp(k)*al1
      if(au .gt. 0.0)goto 14021
      b(k,ic)=0.0
      goto 14031
14021 continue
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)
14031 continue
14011 continue
      d=b(k,ic)-bk
      if(abs(d).le.0.0)goto 14001
      dlx=max(dlx,abs(d))
      if(mm(k) .ne. 0)goto 14051
      nin=nin+1
      if(nin .le. nx)goto 14071
      jxx=1
      goto 14002
14071 continue
      mm(k)=nin
      m(nin)=k
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))
14051 continue
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)
      o=o+d*(xb(k)/xs(k))
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)
14001 continue
14002 continue
      if(jxx.gt.0)goto 13992
      d=svr/xm(0)
      if(d .eq. 0.0)goto 14091
      b(0,ic)=b(0,ic)+d
      dlx=max(dlx,abs(d))
      r=r-d*v
      svr=svr-d*xm(0)
14091 continue
      if(dlx.lt.shr)goto 13992
14100 continue
14101 continue
      nlp=nlp+1
      dlx=0.0
14110 do 14111 l=1,nin
      k=m(l)
      jb=ix(k)
      je=ix(k+1)-1
      jn=ix(k+1)-ix(k)
      bk=b(k,ic)
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))
      gk=dot_product(sc(1:jn),x(jb:je))
      gk=(gk-svr*xb(k))/xs(k)
      u=gk+xv(k,ic)*b(k,ic)
      au=abs(u)-vp(k)*al1
      if(au .gt. 0.0)goto 14131
      b(k,ic)=0.0
      goto 14141
14131 continue
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)
14141 continue
14121 continue
      d=b(k,ic)-bk
      if(abs(d).le.0.0)goto 14111
      dlx=max(dlx,abs(d))
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)
      o=o+d*(xb(k)/xs(k))
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)
14111 continue
14112 continue
      d=svr/xm(0)
      if(d .eq. 0.0)goto 14161
      b(0,ic)=b(0,ic)+d
      dlx=max(dlx,abs(d))
      r=r-d*v
      svr=svr-d*xm(0)
14161 continue
      if(dlx.lt.shr)goto 14102
      goto 14101
14102 continue
      goto 13991
13992 continue
      if(jxx.gt.0)goto 13902
      if(abs(b(0,ic)-bs(0,ic)).gt.shr) ixx=1
      if(ixx .ne. 0)goto 14181
14190 do 14191 j=1,nin
      if(abs(b(m(j),ic)-bs(m(j),ic)) .le. shr)goto 14211
      ixx=1
      goto 14192
14211 continue
14191 continue
14192 continue
14181 continue
      sc=b(0,ic)
      b0=0.0
14220 do 14221 j=1,nin
      l=m(j)
      jb=ix(l)
      je=ix(l+1)-1
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)
      b0=b0-b(l,ic)*xb(l)/xs(l)
14221 continue
14222 continue
      sc=min(max(exmn,sc+b0),exmx)
      sxp=sxp-q(:,ic)
      q(:,ic)=min(max(emin*sxp,exp(dble(sc))),emax*sxp)
      sxp=sxp+q(:,ic)
13901 continue
13902 continue
      if(jxx.gt.0)goto 13892
      if(ixx.eq.0)goto 13892
      if(ig.eq.0)goto 13892
      s=-sum(b(0,:))/nc
      b(0,:)=b(0,:)+s
      sc=s
      b0=0.0
14230 do 14231 j=1,nin
      l=m(j)
      if(vp(l) .gt. 0.0)goto 14251
      s=sum(b(l,:))/nc
      goto 14261
14251 continue
      s=elc(parm,nc,b(l,:),is)
14261 continue
14241 continue
      b(l,:)=b(l,:)-s
      jb=ix(l)
      je=ix(l+1)-1
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)
      b0=b0+s*xb(l)/xs(l)
14231 continue
14232 continue
      sc=sc+b0
      sc=exp(sc)
      sxp=sxp*sc
14270 do 14271 ic=1,nc
      q(:,ic)=q(:,ic)*sc
14271 continue
14272 continue
      nit=nit+1
      if(nit .le. maxit)goto 14291
      jerr=-ilm
      return
14291 continue
      goto 13891
13892 continue
      if(jxx .le. 0)goto 14311
      jerr=-10000-ilm
      goto 13812
14311 continue
      devi=0.0
14320 do 14321 ic=1,nc
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)
      a0(ic,ilm)=b(0,ic)
14330 do 14331 i=1,no
      if(y(i,ic).le.0.0)goto 14331
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))
14331 continue
14332 continue
14321 continue
14322 continue
      kin(ilm)=nin
      alm(ilm)=al
      lmu=ilm
      dev(ilm)=(dev1-devi)/dev1
      if(ig.eq.0)goto 13812
      if(ilm.lt.mnl)goto 13811
      if(flmin.ge.1.0)goto 13811
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13812
      if(dev(ilm).gt.devmax)goto 13812
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13812
13811 continue
13812 continue
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc)
      return
      end
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)
      real a0(nc),ca(nx,nc),x(*),f(nc,n)
      integer ia(*),ix(*),jx(*)
14340 do 14341 ic=1,nc
      f(ic,:)=a0(ic)
14341 continue
14342 continue
14350 do 14351 j=1,nin
      k=ia(j)
      kb=ix(k)
      ke=ix(k+1)-1
14360 do 14361 ic=1,nc
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)
14361 continue
14362 continue
14351 continue
14352 continue
      return
      end
      subroutine psort7 (v,a,ii,jj)
c
c     puts into a the permutation vector which sorts v into
c     increasing order. the array v is not modified.
c     only elements from ii to jj are considered.
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements
c
c     this is a modification of cacm algorithm #347 by r. c. singleton,
c     which is a modified hoare quicksort.
c
      dimension a(jj),v(jj),iu(20),il(20)
      integer t,tt
      integer a
      real v
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=a(ij)
      vt=v(t)
      if (v(a(i)).le.vt) go to 30
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      vt=v(t)
 30   l=j
      if (v(a(j)).ge.vt) go to 50
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      vt=v(t)
      if (v(a(i)).le.vt) go to 50
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      vt=v(t)
      go to 50
 40   a(l)=a(k)
      a(k)=tt
 50   l=l-1
      if (v(a(l)).gt.vt) go to 50
      tt=a(l)
      vtt=v(tt)
 60   k=k+1
      if (v(a(k)).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=a(i+1)
      vt=v(t)
      if (v(a(i)).le.vt) go to 100
      k=i
 110  a(k+1)=a(k)
      k=k-1
      if (vt.lt.v(a(k))) go to 110
      a(k+1)=t
      go to 100
      end
