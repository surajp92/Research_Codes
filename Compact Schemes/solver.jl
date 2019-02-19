include("gauss_seidel.jl")

function solver()
    if flag_multigrid == 1
        if flag_solver == 1
            jacobi_solver()
        elseif flag_solver == 2
            gauss_seidel()
        elseif flag_solver == 3
            steepest_descent()
        elseif flag_solver == 4
            conjugate_gradient()
        else
            biconjugate_gradient()
        end
    else
        multigrid_solver()
    end
end
