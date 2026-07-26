using LinearAlgebra: I

function skew_matrix(vector)
    return [0.0 -vector[3] vector[2];
            vector[3] 0.0 -vector[1];
            -vector[2] vector[1] 0.0]
end

function reconstruct_centerline(node_coordinates, state, flexibility_matrix,
                                initial_curvature)
    state_matrix = reshape(state, 12, :)
    order = sortperm(node_coordinates; alg = Base.Sort.MergeSort)
    x = node_coordinates[order]
    resultants = state_matrix[7:12, order]
    strain_curvature = flexibility_matrix * resultants
    strain = strain_curvature[1:3, :]
    curvature = strain_curvature[4:6, :]

    centerline = zeros(3, length(x))
    rotation = Matrix{Float64}(I, 3, 3)
    e1 = [1.0, 0.0, 0.0]
    for index in 2:length(x)
        step = x[index] - x[index - 1]
        step == 0 && continue
        average_strain = 0.5 *
                         (strain[:, index - 1] + strain[:, index])
        average_curvature = initial_curvature +
                            0.5 *
                            (curvature[:, index - 1] + curvature[:, index])
        rotation_half = rotation *
                        exp(0.5 * step * skew_matrix(average_curvature))
        centerline[:, index] = centerline[:, index - 1] +
                               step * rotation_half *
                               (e1 + average_strain)
        rotation = rotation * exp(step * skew_matrix(average_curvature))
    end

    return x, centerline
end

function root_rotation_angle(time, terminal_angular_speed, ramp_duration)
    if time <= ramp_duration
        return 0.5 * terminal_angular_speed * time^2 / ramp_duration
    else
        return terminal_angular_speed * (time - 0.5 * ramp_duration)
    end
end

function rotation_about_x3(angle)
    cosine = cos(angle)
    sine = sin(angle)
    return [cosine -sine 0.0; sine cosine 0.0; 0.0 0.0 1.0]
end
