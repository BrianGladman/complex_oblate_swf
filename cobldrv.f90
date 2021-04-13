program cobldrv
    use param
    use complex_oblate_swf

    real(knd)       x, ca, arg1, darg, step1, step2, step3, xneu
    complex(knd)    c
    character       chr

    integer         i, im, j, mmin, minc, mnum, m, l, lnum, ioprad, iopang, &
                    iopnorm, ioparg, narg, kind, kindd, kindq

    real (knd), dimension(:), allocatable ::        eta
    complex(knd), dimension(:), allocatable ::      r1c, r1dc, r2c, r2dc
    integer, dimension(:), allocatable ::           ir1e, ir1de, ir2e, ir2de, naccr
    complex(knd), dimension (:,:), allocatable ::   s1c, s1dc
    integer, dimension(:,:), allocatable ::         is1e, is1de, naccs, naccds

    kindd =  8
    kindq = 16

    open(1, file='coblfcn.dat')
    open(20, file='fort.20')
    open(30, file='fort.30')
    open(40, file='fort.20')
    open(50, file='fort.30')
    open(60, file='fort.20')

    read(1,*) mmin, minc, mnum, lnum
        read(1,*) ioprad, iopang, iopnorm
        read(1,*) c, x
        if(iopang.ne.0) read(1,*) ioparg, arg1, darg, narg


    allocate (eta(narg), r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate (ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    allocate (s1c(lnum, narg), s1dc(lnum, narg))
    allocate (is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg), naccds(lnum, narg))

    if(knd.eq.kindd) minacc = 8
    if(knd == kindq .and. aimag(c) <= 20.0e0_knd) minacc = 15
    if(knd == kindq .and. aimag(c) >  20.0e0_knd) minacc =  8

    ca=abs(c)
    ndec=precision(ca)
    nex=range(ca) - 1

    nbp = int(2 * (real(c) + aimag(c)) / 3.1416)
    if(lnum < nbp .and. (2 * (lnum / 2) /= lnum)) lnum = lnum + 1
    maxe = max(50, nbp + 30) + 10
    maxe2 = maxe + maxe
    if(maxe2 < lnum) maxe2 = lnum
    maxm = mmin + minc * (mnum - 1)
    maxlp = lnum + maxm + 1
    maxint = 2 * (nbp + 33) + 306
    maxj = lnum + 3 * ndec + int(ca) + 5 + maxm
    maxp = max(lnum + 3 * ndec + int(ca) + 5, maxlp + 5)
    maxn = maxj
    maxpdr = 2 * int(ca) + 4 * ndec + int(100 * x) + 8
    neta = 30
    step1 = 0.1e0_knd
    step2 = 0.1e0_knd
    step3 = 0.8e0_knd
    nstep1 = 1
    nstep2 = 1
    nstep3 = 3
    ngau = 100 * (2 + int(ca / 500.0e0_knd))
    xneu = 0.3e0_knd

    if(ioprad == 2) then
        if(x < 0.2e0_knd) then
            step1 = x / 4.0e0_knd
            step2 = 0.075e0_knd
            step3 = 1.0e0_knd - step1 - step2
            nstep1 = 1
            nstep2 = 1
            nstep3 = 4
            ngau=100 * (2 + int(ca / 500.0e0_knd))
        end if

        if(x.lt.0.1e0_knd) then
            step1 = x /4.0e0_knd
            step2 = 0.025e0_knd
            step3 = 1.0e0_knd - step1 - step2
            nstep1 = 1
            nstep2 = 1
            nstep3 = 4
            ngau = 200
            if(ca > 500.0e0_knd) ngau = 200 * (2 + int(ca /1000.0e0_knd))
        end if

        if(x < 0.05e0_knd) then
            step1 = x / 15.0e0_knd
            step2 = 0.002e0_knd
            step3 = 1.0e0_knd - step1 - step2
            nstep1 = 1
            nstep2 = 1
            nstep3 = 2
            ngau = 300
            if(ca >  500.0e0_knd) ngau =  500
            if(ca > 1000.0e0_knd) ngau =  800
            if(ca > 1500.0e0_knd) ngau = 1000
            if(ca > 2000.0e0_knd) ngau = 1200
            if(ca > 2500.0e0_knd) ngau = 1500
            if(ca > 3000.0e0_knd) ngau = 1700
            if(ca > 3500.0e0_knd) ngau = 1900
            if(ca > 4000.0e0_knd) ngau = 2200
            if(ca > 4500.0e0_knd) ngau = 2500
            if(ca > 5000.0e0_knd) ngau = 2500 + 300 * int((ca - 4500.0e0_knd) / 500.0e0_knd)
        end if

        if(x < 0.01e0_knd) then
            step1 = x / 15.0e0_knd
            step2 = 0.002e0_knd
            step3 = 1.0e0_knd - step1 - step2
            nstep1 = 1
            nstep2 = 1
            nstep3 = 2
            ngau = 300
            if(ca >  500.0e0_knd) ngau =  600
            if(ca > 1000.0e0_knd) ngau =  800
            if(ca > 1500.0e0_knd) ngau = 1000
            if(ca > 2000.0e0_knd) ngau =  300 * int(ca /500.0e0_knd)
        end if

        if(x < 0.001e0_knd) then
            step1 = x / 15.0e0_knd
            step2 = 0.002e0_knd
            step3 = 1.0e0_knd - step1 - step2
            nstep1 = 3
            nstep2 = 1
            nstep3 = 2
            ngau = 600
            if(ca >  300.0e0_knd) ngau =  800
            if(ca > 1000.0e0_knd) ngau =  900
            if(ca > 1500.0e0_knd) ngau = 1000
            if(ca > 2000.0e0_knd) ngau = 300 * int(ca / 500.0e0_knd)
            if(aimag(c) > 5.0e0_knd) ngau = ngau + (ngau / 5) * min(5, int(aimag(c)) - 5)
            if(knd.eq.kindq.and.x.lt.0.00001e0_knd) ngau=2*ngau
            if(knd.eq.kindq.and.x.lt.0.000001e0_knd) ngau=2*ngau
            if(knd.eq.kindd.and.aimag(c).gt.4.0e0_knd.and.x.lt. 0.00001e0_knd) ngau=2*ngau
            if(knd.eq.kindd.and.aimag(c).gt.4.0e0_knd.and.x.lt. 0.000001e0_knd) ngau=2*ngau
        end if

        xneu = 0.3e0_knd
        if(ca > 100.0e0_knd) xneu = 0.04e0_knd
        if(ca > 600.0e0_knd) xneu = 0.03e0_knd
        if(ca > 800.0e0_knd) xneu = 0.01e0_knd
        if(aimag(c) > 50.0e0_knd) xneu = 0.01e0_knd

        if(x >= 0.01e0_knd) then
            if(knd == kindd) then
                if(x >= 0.01e0_knd) maxn = 2 * int(25 / (x * x) + 300 / x + 3 * ca + 1250 * knd) + 5
                if(x >= 0.1e0_knd)  maxn = 2 * int((lnum + ca / 5 + 0.5e0_knd * maxm + 200) * 1.4e0_knd / x) + 5
                if(x >= 0.5e0_knd)  maxn = 2 * int((lnum + ca / 5 + 0.5e0_knd * maxm + 300) / x) + 5
                if(x >= 1.0e0_knd)  maxn = 2 * int(lnum + ca / 5 + 0.5e0_knd * maxm + 300) + 5
            end if
            if(knd.eq.kindq) then
                if(x >= 0.01e0_knd) maxn = 2 * int(25 / (x * x) + 400 / x + 3 * ca + 1250 * knd) + 5
                if(x >= 0.1e0_knd)  maxn = 2 * int((lnum + ca / 5 + 0.5e0_knd * maxm+ 350) * 1.4e0_knd / x) + 5
                if(x >= 0.5e0_knd)  maxn = 2 * int((lnum + ca / 5 + 0.5e0_knd * maxm + 400) / x) + 5
                if(x >= 1.0e0_knd)  maxn = 2 * int(lnum + ca / 5 + 0.5e0_knd * maxm + 400) + 5
                end if
            maxn = maxn + maxm
        end if
    end if

    maxp = max(maxn, maxp, maxpdr)
    maxq = lnum + 3 * ndec + int(ca) + maxm + maxm + 4
    if(knd == kindd .and. aimag(c) < 10.0e0_knd .and. real(c) <= 60.0e0_knd &
        .and. real(c) >= 10.0e0_knd .and. mmin <= 40 .and. x <= 0.99e0_knd &
        .and. x > 0.1e0_knd) maxq = max(maxq, 250 - int(50 * x) + 4 + maxm + maxm)
    maxdr = maxpdr / 2 + 1
    maxd = maxn / 2 + 1
    maxmp = maxm + maxm + 5
    maxt = 1
    jnebmax = 30
    jnenmax = 10
    if(x < 0.05e0_knd) jnenmax = 1
    if(iopang /= 0) maxt = narg

    if (iopang /= 0) then
        do j = 1, narg
            eta(j) = arg1 + (j - 1) * darg
        end do
    end if

    do im = 1, mnum
        m = mmin + (im - 1) * minc

        call coblfcn(c, m, lnum, ioprad, x, iopang, iopnorm, narg, eta, &
                     r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                     s1c, is1e, s1dc, is1de, naccs, naccds)

        if (ioprad /= 0) then

            if(knd == kindd) write(20, 10) x, c, m
10          format(1x,e23.14,e23.14,e23.14,i5)
            if(knd == kindq) write(20,20) x, c, m
20          format(1x,e39.30,e39.30,e39.30,i5)

            do i = 1, lnum
                l = m + i - 1
                chr = 'w'
                if (naccr(i) < 0) chr = 'e'
                write(20, 30) l, r1c(i), ir1e(i), r1dc(i), ir1de(i), r2c(i), ir2e(i), r2dc(i), ir2de(i), abs(naccr(i)), chr
30              format(1x,i5,2x,2(f17.14,1x,f17.14,i6,2x),/,8x,2(f17.14,1x,f17.14,i6,2x),i2,a)
            end do
        end if

        if (iopang /= 0) then

            if(knd == kindd) write(30, 40) c, m
40          format(1x,'c = ',e23.14,e23.14,'; m = ',i5)
            if(knd == kindq) write(30, 50) c, m
50          format(1x,'c = ',e39.30,e39.30,'; m = ',i5)

            do i = 1, lnum
                l = m + i - 1

                write(30, "(1x,i6)") l

                do j = 1, narg
                    arg = arg1 + (j - 1) * darg

                    if(ioparg.eq.1.and.iopang.eq.1) write(30, 60) eta(j), s1c(i,j), is1e(i,j), naccs(i,j)
                    if(ioparg.eq.1.and.iopang.eq.2) write(30, 70) eta(j), s1c(i,j), is1e(i,j), s1dc(i,j), is1de(i,j), naccs(i,j), naccds(i,j)

60                  format(1x,f19.14,2x,f17.14,1x,f17.14,2x,i5,2x,', ',i2)
70                  format(1x,f19.14,2x,f17.14,1x,f17.14,2x,i5,2x,f17.14,1x,f17.14,2x,i5,2x,i2,', ',i2)

                end do
            end do
        end if
    end do

    deallocate (is1e, is1de, naccs, naccds)
    deallocate (s1c, s1dc)
    deallocate (ir1e, ir1de, ir2e, ir2de, naccr)
    deallocate (eta, r1c, r1dc, r2c, r2dc)
    close(60)
    close(50)
    close(40)
    close(30)
    close(20)
    close(1)

end program cobldrv
