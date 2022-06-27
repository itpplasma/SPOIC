using HDF5
using Polynomials
using OffsetArray

function get_zernike(r, lrad, mpol, zernike)
  rm = 1.0
  rm1 = 0.0
  zernike[:,:,:] = 0.0
  for m in range(0, mpol)
    if (lrad >= m)
      zernike[m,m,0:1] = [rm, real(m)*rm1]
    end

    if (lrad >= m+2) then
      zernike[m+2,m,0] = real(m+2)*rm*r^2 - real(m+1)*rm
      zernike[m+2,m,1] = real((m+2)^2)*rm*r - real((m+1)*m)*rm1
    end

    for n in range(m+4, lrad, step=2)
      factor1 = real(n)/real(n^2 - m^2)
      factor2 = real(4 * (n-1))
      factor3 = real((n-2+m)^2)/real(n-2) + real((n-m)^2)/real(n)
      factor4 = real((n-2)^2-m^2) / real(n-2)

      zernike[n, m, 0] = factor1 * ((factor2*r^2 - factor3)*zernike[n-2,m,0] - factor4*zernike[n-4,m,0])
      zernike[n, m, 1] = factor1 * (two*factor2*r*zernike[n-2,m,0] + (factor2*r^2 - factor3)*zernike[n-2,m,1] - factor4*zernike[n-4,m,1])
    end

    rm1 = rm
    rm = rm * r

  end

  for n in range(2, lrad, step=2)
    zernike[n,0,0] = zernike[n,0,0] - (-1)^(n/2)
  end

  if (mpol >= 1)
    for n in range(3, lrad, step=2)
      zernike[n,1,0] = zernike[n,1,0] - (-1)^((n-1)/2) * real((n+1)/2) * r
      zernike[n,1,1] = zernike[n,1,1] - (-1)^((n-1)/2) * real((n+1)/2)
    end
  end

  for m in range(0, mpol)
    for n in range(m, lrad, step=2)
      zernike[n,m,:] = zernike[n,m,:] / real(n+1)
    end
  end
end

fid = h5open("spec_out.h5", "r")

Ate = fid["vector_potential/Ate"][:,:]
Aze = fid["vector_potential/Aze"][:,:]
Mpol = fid["input/physics/Mpol"]
Lrad = fid["input/physics/Lrad"][:]
Mvol = fid["output/Mvol"][1]
mn = fid["output/mn"][1]
m = fid["output/im"][:]
n = fid["output/in"][:]

close(fid)

zernike = OffsetArray{Float64}(undef, 0:Lrad(1), 0:Mpol, 0:1)

for vvol in range(1,Mvol)
    println("Volume: ", vvol)
    teta = 0.5
    zeta = 0.3
    for ii in range(1, mn)
        mi = m[ii]
        ni = n[ii]
        arg = mi * teta - ni * zeta
        carg = cos(arg)
        sarg = sin(arg)
        if (vvol == 1)
            zernike[ll,mi,0]  # TODO
        else
            T = ChebyshevT(Ate[:])
        end
    end
end
