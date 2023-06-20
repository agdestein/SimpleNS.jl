"""
    process!(processor, stepper)

Process iteration.
"""
function process! end

function process!(logger::Logger, stepper)
    (; n, t) = stepper
    @printf "Iteration %d\tt = %g\n" n t
    logger
end

function process!(o::StateObserver, stepper)
    (; V, p, t) = stepper
    o.state[] = (V, p, t)
end

function process!(writer::VTKWriter, stepper)
    (; setup, V, p, t) = stepper
    (; grid) = setup
    (; xp, yp) = grid

    coords = (xp, yp)

    tformat = replace(string(t), "." => "p")
    vtk_grid("$(writer.dir)/$(writer.filename)_t=$tformat", coords...) do vtk
        vels = get_velocity(V, setup)

        # ParaView prefers 3D vectors. Add zero z-component.
        wp = zeros(size(vels[1]))
        vels = (vels..., wp)
        vtk["velocity"] = vels
        vtk["pressure"] = p
        writer.pvd[t] = vtk
    end

    writer
end
