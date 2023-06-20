"""
    save_vtk(V, p, t, setup, filename = "output/solution")

Save velocity and pressure field to a VTK file.

In the case of a 2D setup, the velocity field is saved as a 3D vector with a
z-component of zero, as this seems to be preferred by ParaView.
"""
function save_vtk(V, p, setup, filename = "output/solution")
    parts = split(filename, "/")
    path = join(parts[1:end-1], "/")
    isdir(path) || mkpath(path)
    (; xp, yp) = setup.grid
    vtk_grid(filename, xp, yp) do vtk
        vels = get_velocity(V, setup)

        # ParaView prefers 3D vectors. Add zero z-component.
        wp = zeros(size(vels[1]))
        vels = (vels..., wp)
        vtk["velocity"] = vels
        vtk["pressure"] = p
    end
end
