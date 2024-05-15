## Julia set

$z$ is repeatedly updated using:  $z \leftarrow z^2 + c$  where c is another complex number that gives a specific Julia set. After numerous iterations, if the magnitude of z is less than 2 we say that pixel is in the Julia set. Performing this calculation for a whole grid of pixels gives a fractal image. 



```
make TARGET=JuliaSet.cpp
make run
```



```
convert -resize 100% -delay 20 -loop 0 *.png JuliaSet.gif
```



## Mandelbrot set

The Mandelbrot set is the set of all $c$ for which the iteration $z \leftarrow z^2 + c$, starting from $z = 0$, does not diverge to infinity. Julia sets are either connected (one piece) or a dust of infinitely many points. The Mandelbrot set is those c for which the Julia set is connected.


```
make TARGET=MandlebrotSet.cpp
make run
```

