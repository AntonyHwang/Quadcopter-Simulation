%simulation times, in seconds.
start_time = 0;
end_time = 10;
dt = 0.005;
times = start_time:dt:end_time;

%number of points in the simualation.
N = numel(times);

%initial simulation state.
x = [0; 0; 10];
xdot = zeros(3, 1);
theta = zeros(3, 1);

%simulate some disturbance in the angular velocity.
%the magnitude of the deviation is in radians/second.
deviation = 100;
thatadot = deg2rad(2 * deviation * rand(3, 1) - deviation);

%step through the simulation, updating the state.
for t = times
    %take input form our sontroller.
    i = input(t);
    
    omega = thetadot2omega(thetadot, theta);
    
    %compute linear and angular accelerations.
    a = acceleration(i, theta, xdot, m, g, k, kd);
    omegadot = angular_acceleration(i, omega, I, L, b, k);
    
    omega = omega + dt * omegadot;
    thetadot = omega2thetadot(omega, theta);
    theta = theta + dt * thetadot;
    xdot = xdot + dt * a;
    x = x + dt * xdot;
end

%compute thrust given current inputs and thrust coefficient.
function T = thrust(inputs, k)
    %inputs are values for w_i^2
    T = [0; 0; k * sum(inputs)];
end

%compute torques, given current inputs, length, drag coefficient, and
%thrust coefficient
function tau = torques(inputs, L, b, k)
    %inputs are values for w_i^2
    tau = [L * k * (inputs(1) - inputs(3));
           L * k * (inputs(2) - inputs(4));
           b * (inputs(1) - inputs(2) + inputs(3) - inputs(4))];
end

function a = acceleration(inputs, angles, xdot, m, g, k, kd)
    gravity = [0; 0; -g];
    R = rotation(angles);
    T = R * thrust(inputs, k);
    Fd = -kd * xdot;
    a = gravity + 1 / m * T + Fd;
end

function omegadot = angular_acceleration(inputs, omega, I, L, b, k)
    tau = torques(inputs, L, b, k);
    omegadot = inv(I) * (tau - cross(omega, I * omega));
end

%compute system inputs and updated state.
function [input, state] = pd_controller(state, thetadot)
    %controller gains, tuned by hand and intuition.
    Kd = 4;
    Kp = 3;
    
    %initialize the integral if necessary.
    if isfield(state, 'integral')
        params.integral = zeros(3, 1);
    end
    
    %compute total thrust.
    total = state.m * state.g / state.k / (cos(state.integral(1)) * cos(state.integral))
    
    %compute errors.
    e = Kd * thetadot + Kp * params.integral;
    
    %solve for the inputs, gamma_i.
    input = error2inputs(params, accels, total);
    
    %update the state.
    params.integral = params.integral + params.dt .* thetadot;
end

%compute system inputs and updated state.
function [input, state] = pid_controller(state, thetadot)
    %controller gains, tuned by hand and intuition.
    Kd = 4;
    Kp = 3;
    Ki = 5.5;
    
    %initialize the integral if necessary.
    if isfield(state, 'integral')
        params.integral = zeros(3, 1);
        params.integral2 = zeros(3, 1);
    end
    
    %prevent wind-up.
    if max(abs(params.integral2)) > 0.01
        params.integral2(:) = 0;
    end
    
    %compute total thrust.
    total = state.m * state.g / state.k / (cos(state.integral(1)) * cos(state.integral))
    
    %compute errors.
    e = Kd * thetadot + Kp * params.integral - Ki * params.integral2;
    
    %solve for the inputs, gamma_i.
    input = error2inputs(params, accels, total);
    
    %update the state.
    params.integral = params.integral + params.dt .* thetadot;
    params.integral2 = params.integral2 + params.dt .* thetadot;
end

function J = cost(theta)
    %create a controller using the given gains.
    control = controller('pid', theta(1), theta(2), theta(3));
    
    %perform simulation.
    data = simulate(control);
    
    %compute the integral.
    t0 = 0;
    tf = 1;
    J = 1 / (tf - t0) * sum(data.theta(data.t >= t0 & data.t <= tf) .^ 2) * data.dt;
end
