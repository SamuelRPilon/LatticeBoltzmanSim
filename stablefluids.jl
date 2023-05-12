## SAMUEL R PILON
#JULIALANG - FLuid Simulation using lattice boltzmann method
# For personal use of SAMUEL PILON only. Not to be distributed

"""
__CODE DESCRIPTION__

The Lattice Boltzmann Simulation can be solved using FFT to obtain ultra-fast simulations. 
The simulation involves solving the momentum and incompressibility equations. The momentum 
equation includes the velocity, pressure, forcing, kinematic viscosity, time, and the nabla 
and Laplace operators.

The solution strategy involves starting with zero velocity everywhere and then adding forces 
to the system. The forces are convected by self-advection, where the value at the current 
location is set to be the value at the position backtraced on the streamline. This step is 
unconditionally stable.

Next, the system is diffused and projected in the Fourier domain using a series of steps.
First, the forward transformation is applied to the convected velocity in the Fourier domain.
Then, the velocity is diffused using low-pass filtering, which involves convolution in the Fourier
domain. The (pseudo-) pressure is computed in the Fourier domain by evaluating the divergence,
and the velocities are corrected to ensure incompressibility. Finally, the inverse transformation 
is applied to bring the velocities back into the spatial domain.

This process is repeated to obtain the simulation results. The spatial frequencies (wavenumbers) 
are represented by k = [kₓ, k_y], and the Fourier transformation implicitly prescribes periodicity.

%%%%%%%%%%%%%%%%%%%%%%%%%
-PACKAGES NEEDED 

* FFTW  
* Plots 
* ProgressMeter 
* Interpolations
* LinearAlgebra 

%%%%%%%%%%%%%%%%%%%%%%%%%

"""
# Imports Libraries
using FFTW  # Fast Fourier Transform
using Plots # Self Explainatory
using ProgressMeter # Adds a Progress bar to display in Terminal
using Interpolations # Makes Interpolation Easier
using LinearAlgebra # Some math Library

# Intitial Conditions
numpoints = 250
mu = 0.0001
tstep = 0.01
tspan = 300

# the Fucntion utilizes all the other functions and runs in the main fuctions to compile
function main()
    elength = 1.0 / (numpoints - 1)
    xint = 0.0:elength:1.0
    y_interval = 0.0:elength:1.0

    # Meshes the grid for the calculation to occurr
    coordinates_x = [x for x in xint, y in y_interval]
    coordinates_y = [y for x in xint, y in y_interval]

    wavenumbers_1d = fftfreq(numpoints) .* numpoints

    wavenumbers_x = [k_x for k_x in wavenumbers_1d, k_y in wavenumbers_1d]
    wavenumbers_y = [k_y for k_x in wavenumbers_1d, k_y in wavenumbers_1d]
    wavenumbers_norm = [norm([k_x, k_y]) for k_x in wavenumbers_1d, k_y in wavenumbers_1d]

    decay = exp.(- tstep .* mu .* wavenumbers_norm.^2)

    wavenumbers_norm[iszero.(wavenumbers_norm)] .= 1.0
    normalized_wavenumbers_x = wavenumbers_x ./ wavenumbers_norm
    normalized_wavenumbers_y = wavenumbers_y ./ wavenumbers_norm

    # Define the forces that act on the fluid for the inital motions to occurr
    force_x = 100.0 .* (
        exp.( - 1.0 / (2 * 0.005) * ((coordinates_x .- 0.2).^2 + (coordinates_y .- 0.45).^2))
        -
        exp.( - 1.0 / (2 * 0.005) * ((coordinates_x .- 0.8).^2 + (coordinates_y .- 0.55).^2))
    )

    # Makes all the data frames or i think thats what they are called. 
    backtraced_coordinates_x = zero(coordinates_x)
    backtraced_coordinates_y = zero(coordinates_y)

    velocity_x = zero(coordinates_x)
    velocity_y = zero(coordinates_y)

    velocity_x_prev = zero(velocity_x)
    velocity_y_prev = zero(velocity_y)

    velocity_x_fft = zero(velocity_x)
    velocity_y_fft = zero(velocity_y)
    pressure_fft = zero(coordinates_x)

    # change the dark if you want a dark theme, the light theme mathes better with my website. 
    theme(:default)

    @showprogress "Timestepping ..." for iter in 1:tspan

        # Apply the forces that were initilized and defined 
        time_current = (iter - 1) * tstep
        pre_factor = max(1 - time_current, 0.0)
        velocity_x_prev += tstep * pre_factor * force_x

        # The fluid knows where it is at all times. It knows this because it knows where 
        #it isn't, by subtracting where it is, from where it isn't, or where it isn't, from 
        #where it is, whichever is greater, it obtains a difference, or deviation. 
        #The guidance sub-system uses deviations to generate corrective commands to drive 
        #the fluid from a position where it is, to a position where it isn't, and arriving 
        #at a position where it wasn't, it now is. Consequently, the position where it is, is 
        #now the position that it wasn't, and it follows that the position where it was, is now 
        #the position that it isn't
        backtrace!(backtraced_coordinates_x, coordinates_x, velocity_x_prev)
        backtrace!(backtraced_coordinates_y, coordinates_y, velocity_y_prev)
        interpolate_positions!(
            velocity_x,
            velocity_x_prev,
            xint,
            y_interval,
            backtraced_coordinates_x,
            backtraced_coordinates_y,
        )
        interpolate_positions!(
            velocity_y,
            velocity_y_prev,
            xint,
            y_interval,
            backtraced_coordinates_x,
            backtraced_coordinates_y,
        )

        # Transforms into a domain in which the fourier series can be calulcated 
        velocity_x_fft = fft(velocity_x)
        velocity_y_fft = fft(velocity_y)

        # stack overflow said that this is needed to get rid of the weird aliasing that was happening
        velocity_x_fft .*= decay
        velocity_y_fft .*= decay

        # Compute Pseudo-Pressure by Divergence in Fourier Domain (they also said to in do this)
        pressure_fft = (
            velocity_x_fft .* normalized_wavenumbers_x
            +
            velocity_y_fft .* normalized_wavenumbers_y
        )

        # they are incompressible (307 taught me this)
        velocity_x_fft -= pressure_fft .* normalized_wavenumbers_x
        velocity_y_fft -= pressure_fft .* normalized_wavenumbers_y

        # transform out of fourier domain into a 2d spatial domain
        velocity_x = real(ifft(velocity_x_fft))
        velocity_y = real(ifft(velocity_y_fft))

        # spans the animation
        velocity_x_prev = velocity_x
        velocity_y_prev = velocity_y

        # Holds the plot and spans the animation for the tspan
        d_u__d_y = diff(velocity_x, dims=2)[2:end, :]
        d_v__d_x = diff(velocity_y, dims=1)[:, 2:end]
        curl = d_u__d_y - d_v__d_x
        display(heatmap(xint, y_interval, curl', c=:plasma, aspect_ratio=:equal, size=(1000, 1000),title=("Lattice Boltzmann Simulation"), titlefont=font(18, :center, :bottom)))

    end
end

# The function first constructs a LinearInterpolation object using the field 
# data and the grid point coordinates x and y. This creates a 
# function that can interpolate values of the scalar field at any arbitrary position 
# within the bounds of the grid.

function interpolate_positions!(
    interp_field, field, x, y, qpx, qpy,
)
    interpolator = LinearInterpolation(
        (x, y),
        field,
    )
    interp_field[:] = interpolator.(qpx, qpy)
end

    # The function uses the Euler method to step backwards in time by a
    # constant time step length (tstep) and periodically clamps 
    # the resulting positions to lie within the domain [0.0, 1.0].
    
    function backtrace!(
        back_pos,
        ogpos,
        direction,
    )
        # Euler Step backwards in time and periodically clamp into [0.0, 1.0]
        # The mod1() function in Julia is used to perform periodic wrapping of the 
        # positions, ensuring that they remain within the domain
        back_pos[:] = mod1.(
            ogpos - tstep * direction,
            1.0,
        )
    end

main()