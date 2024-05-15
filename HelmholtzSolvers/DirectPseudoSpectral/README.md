## One-Dimensional Helmholtz Equation in the Cartesian Coordinate System

We solve the following equation using the Finite Difference Method(FDM)

$$\frac{d^2 u}{dx^2} - \lambda u = f,\quad \lambda \geq 0$$

with periodic boundary condition.

To build and run the code, invoke the following commands
```
make TARGET=TestAllPeriodic1D.cpp
make run
```

## Two-Dimensional Helmholtz Equation in the Cartesian Coordinate System

We solve the following equation using the Finite Difference Method(FDM)

$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} - \lambda u = f,\quad \lambda \geq 0$$

with the periodic boundary condition along the x- and y-directions.

To build and run the code, invoke the following commands

```
make TARGET=TestAllNonPeriodic2D.cpp
make run
```

## Three-Dimensional Helmholtz Equation in the Cartesian Coordinate System

We solve the following equation using the Finite Difference Method(FDM)

$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2} - \lambda u = f,\quad \lambda \geq 0$$

with the periodic boundary condition along the x-, y- and z-directions.

To build and run the code, invoke the following commands

```
make TARGET=TestAllPeriodic3D.cpp
make run
```
