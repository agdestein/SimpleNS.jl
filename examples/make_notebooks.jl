# Convert example scripts to Jupyter notebooks

using Literate

# Scripts to convert
examples = [
    "DecayingTurbulence2D.jl"
    "TaylorGreenVortex2D.jl"
]

# Location of notebooks
output_dir = "notebooks"
ispath(output_dir) || mkpath(output_dir)

# Convert scripts (set `execute = true` to run notebooks upon conversion)
for e ∈ examples
    Literate.notebook(e, output_dir; execute = false)
end
