% 3D LAPLACIAN FUNCTION (FINITE DIFFERENCE)

function lap = laplacian_3d(u)

    % Use periodic-style shifts (equivalent to np.roll)
    lap = -6 * u ...
        + circshift(u, [1, 0, 0]) + circshift(u, [-1, 0, 0]) ...
        + circshift(u, [0, 1, 0]) + circshift(u, [0, -1, 0]) ...
        + circshift(u, [0, 0, 1]) + circshift(u, [0, 0, -1]);

end