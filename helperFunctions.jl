using LinearAlgebra
using SparseArrays



function compLocalMin(time)
    derivative = diff(diff(time))                       # compute the first derivative of timstep sizes (FD)
    
    dt = Float64[]
    

    for i in 1:length(derivative)-1                     # Check where the derivative changes sign
            
        if derivative[i] <= 0 && derivative[i+1] >= 0   # Upcrossing: local minimum
            
            push!(dt, diff(time)[i+1])                  # save the value of stepsize(NOT derivative)
            
        end
    end

    return  dt
end

function nLV!(A, Ny, η_eff,h)
    
    A[1,1]      = 1
    A[end,end]  = 1
    for i in 2:Ny-1
        A[i,i+1] = η_eff[i]/2h^2
        A[i,i] = -(η_eff[i-1] + η_eff[i])/2h^2
        A[i,i-1] = η_eff[i-1]/2h^2
    end
    

end



function D2nd(n, Lx)
    x = range(0,Lx,n+1) 
    xc = (x[2:end] + x[1:end-1]) ./ 2
    h = abs(xc[2]-xc[1])

    du = ones(n-1)
    dl = ones(n-1) 
    d = ones(n) .*-2
    D = Tridiagonal(dl,d,du)

    D = D .* (1/(h^2))

    return xc, D, h

end