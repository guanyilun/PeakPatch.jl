# Time Healpix.jl map2alm / alm2cl at Nside 4096 on a released Websky map.
using Healpix, Printf
m = Healpix.readMapFromFITS("/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/ksz.fits", 1, Float64)
@info "read" nside=m.resolution.nside threads=Threads.nthreads()
for lmax in (2048, 8192)
    t = @elapsed alm = Healpix.map2alm(m; lmax=lmax, mmax=lmax, niter=0)
    t2 = @elapsed cl = Healpix.alm2cl(alm)
    @printf("lmax=%d: map2alm %.1f s, alm2cl %.2f s; D(3000)=%.3f uK2\n", lmax, t, t2,
            lmax >= 3000 ? 3000*3001/2π*cl[3001] : NaN)
end
