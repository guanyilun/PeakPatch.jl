#!/usr/bin/env julia
# Measure spectral moments sigma0,1,2 of Julia's Rf=2.068-smoothed field vs theory
# (1.783, 1.280, 3.703). If sigma2 is low, Julia's field lacks near-Nyquist power
# -> fewer grid-scale peaks at the smallest filter (the resolution-floor halos).
using PeakPatch, FFTW, Printf, Statistics
import PeakPatch.PowerSpectrum: load_pk
import PeakPatch.Filters: smooth_field, tophat_window

n=256; box=320.8333; seed=13579; Rf=2.068
pk=load_pk(joinpath(@__DIR__,"data","pk_websky.dat"))
delta=generate_grf(n,pk,box,seed; fortran_compat=true)
delta_k=rfft(delta)

# sigma0: smoothed field RMS
ds=smooth_field(delta_k,n,box,Rf,1; fortran_compat=true)
s0=sqrt(mean(ds.^2))

# sigma1^2 = <|grad d_s|^2> = int k^2 P W^2 ; sigma2^2=<(lap d_s)^2>=int k^4 P W^2
dk=2pi/box
kx=dk.*Float64.(FFTW.rfftfreq(n,n)); ky=dk.*Float64.(FFTW.fftfreq(n,n)); kz=dk.*Float64.(FFTW.fftfreq(n,n))
sm_k=copy(delta_k)
for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
    kk=sqrt(kx[ix]^2+ky[iy]^2+kz[iz]^2)
    sm_k[ix,iy,iz]*= kk==0 ? 0.0 : tophat_window(kk*Rf)
end
# grad^2 and lap^2 via k-space
s1sq=0.0; s2sq=0.0; norm=Float64(n)^6  # for power normalization of rfft
# Easier: build fields. lap = irfft(-k^2 * sm_k)
lapk=similar(sm_k);
for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
    kk2=kx[ix]^2+ky[iy]^2+kz[iz]^2
    lapk[ix,iy,iz]=sm_k[ix,iy,iz]*(-kk2)
end
lap=irfft(lapk,n)
s2=sqrt(mean(lap.^2))
# grad: |grad|^2 sum of three derivative fields
gsq=zeros(Float64,n,n,n)
for (dimk,kv) in enumerate((kx,ky,kz))
    gk=similar(sm_k)
    for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
        kcomp = dimk==1 ? kx[ix] : dimk==2 ? ky[iy] : kz[iz]
        gk[ix,iy,iz]=sm_k[ix,iy,iz]*im*kcomp
    end
    g=irfft(gk,n); gsq.+=g.^2
end
s1=sqrt(mean(gsq))
@printf("Julia smoothed-field moments at Rf=%.3f:\n",Rf)
@printf("  sigma0 = %.4f   (theory 1.783)\n",s0)
@printf("  sigma1 = %.4f   (theory 1.280)\n",s1)
@printf("  sigma2 = %.4f   (theory 3.703)\n",s2)
@printf("  gamma=s1^2/(s0 s2)=%.4f  R*=sqrt3 s1/s2=%.4f Mpc/h (theory R*=0.599)\n",
        s1^2/(s0*s2), sqrt(3)*s1/s2)
