# Run: julia --project=validation -t 8 validation/quantify_2lpt_nyquist.jl   (see NOTES_2LPT_NYQUIST_2026-09-24.md)
# Quantify the production GPU tile-local 2LPT Nyquist mismatch vs the exact determinant.
using PeakPatch, FFTW, Statistics, Printf
FFTW.set_num_threads(Threads.nthreads())
const LPT = PeakPatch.LPT
pk = PeakPatch.PowerSpectrum.load_pk("validation/websky_6144/data/pk_websky.dat")
for n in (128, 256)
    L = n * 0.85221                      # production cell size (Mpc/h)
    δ = Float64.(PeakPatch.RandomField.generate_grf(n, pk, L, 12345))
    δk = rfft(δ)
    ψ1 = LPT.displacements_1lpt(δk, n, L)
    ψ2 = LPT.displacements_2lpt(δk, n, L)                       # exact (serial convention)
    dk = 2π / L; kx = rfftfreq(n, n * dk); ky = fftfreq(n, n * dk); nk = n ÷ 2 + 1; nyq = n ÷ 2
    bad(ix, iy, iz) = ix == nk || iy == nyq + 1 || iz == nyq + 1
    # production convention: φ_ij zero k=0 AND Nyquist, trace term = raw δ²
    src = δ .^ 2 .* 0.5
    for (i, j, c) in ((1,1,.5),(2,2,.5),(3,3,.5),(1,2,1.),(1,3,1.),(2,3,1.))
        φk = similar(δk)
        for iz in 1:n, iy in 1:n, ix in 1:nk
            k = (kx[ix], ky[iy], ky[iz]); k2 = sum(abs2, k)
            φk[ix,iy,iz] = (k2 == 0 || bad(ix,iy,iz)) ? 0 : -k[i]*k[j]/k2 * δk[ix,iy,iz]
        end
        src .-= irfft(φk, n) .^ 2 .* c
    end
    sk = rfft(src)
    ψ2m = map(1:3) do d
        o = similar(sk)
        for iz in 1:n, iy in 1:n, ix in 1:nk
            k = (kx[ix], ky[iy], ky[iz]); k2 = sum(abs2, k)
            o[ix,iy,iz] = (k2 == 0 || bad(ix,iy,iz)) ? 0 : -im*k[d]/k2 * sk[ix,iy,iz]
        end
        irfft(o, n)
    end
    smooth(f, R) = R == 0 ? f : irfft(rfft(f) .* [exp(-0.5*R^2*(kx[ix]^2+ky[iy]^2+ky[iz]^2)) for ix in 1:nk, iy in 1:n, iz in 1:n], n)
    @printf("n=%d (L=%.0f Mpc/h, Nyquist-plane mode fraction %.2f%%)\n", n, L, 100*(3n^2)/(n^3/2))
    for R in (0.0, 2.0, 5.0)
        e = sqrt(sum(mean(abs2, smooth(ψ2m[d] .- ψ2[d], R)) for d in 1:3))
        a2 = sqrt(sum(mean(abs2, smooth(ψ2[d], R)) for d in 1:3))
        a1 = sqrt(sum(mean(abs2, smooth(ψ1[d], R)) for d in 1:3))
        @printf("  R=%3.0f Mpc/h: rms Δψ₂/ψ₂ = %.4f ; z=0 total-displacement error (3/7)Δψ₂/ψ₁ = %.2e\n", R, e/a2, (3/7)*e/a1)
    end
end
