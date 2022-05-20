# Robust Quantization for Polycube-Maps

![](figure/rb_teaser.png)

[François Protais](fprotais.github.op), [Maxence Reberol](https://mxncr.github.io/), Nicolas Ray, [Etienne Corman](https://members.loria.fr/ECorman/), Franck Ledoux and [Dmitry Sokolov](https://members.loria.fr/DSokolov/)

# TODO
- Complete this README: detail every binary with inputs and outputs. 
- Splitting potentially bad tets: something similar to [this](https://github.com/fprotais/preprocess_polycube), but more efficient
- Finish incomplete functions
- Improve handling graphite and LinearSolver 

# Future works
- Ovelaps in input polycuboid:
	- Improved voxelisation
	- Deformation that is only injective, without tetgen (add exiting normal constraints)
	- Check that function are actually compatible
- Improve flagging: link with other works or graphcut?
- Custom integer solver to remove MILP solvers



# Use CMake to build the project:
```sh
git clone --recurse-submodules https://github.com/fprotais/robustPolycube &&
cd robustPolycube &&
mkdir build &&
cd build &&
cmake -DCMAKE_BUILD_TYPE=Release .. &&
make -j 
```

# Running the code :

coming...
```sh
./rb_fromscratch ../mesh/S3.mesh 1. hexmesh.mesh 
```
For the supported mesh formats, see [ultimaille](https://github.com/ssloy/ultimaille). 
