# Dinámica Molecular con Termostato y Baróstato de Nosé–Hoover

Este repositorio forma parte de mi Trabajo Fin de Grado titulado **“Simulación de la Entropía en sistemas clásicos”**, presentado en la Universidad de Salamanca. El objetivo general del proyecto es ilustrar el aumento de la entropía en sistemas físicos a través de simulaciones computacionales. En particular, esta sección del trabajo se centra en el estudio de un fluido denso modelado mediante el potencial de Lennard-Jones, bajo condiciones termodinámicas controladas por termostatos y baróstatos de Nosé–Hoover.

## Archivos principales

- `termotato.f95`: Código Fortran que implementa un **termóstato de Nosé–Hoover** para simular un sistema bajo condiciones del ensamble canónico (NVT). Permite estudiar procesos como el enfriamiento lento y la relajación térmica en sistemas de partículas con interacción de tipo Lennard-Jones.

- `barostato.f90`: Código Fortran que extiende el sistema anterior para incluir un **baróstato de Nosé–Hoover**, permitiendo simular en el ensamble isobárico-isotérmico (NPT). Este código ajusta dinámicamente el volumen de la caja de simulación para mantener la presión constante, permitiendo la observación de fenómenos como la cristalización inducida por enfriamiento y compresión.


## Objetivo computacional

Estos programas permiten simular de forma realista la evolución dinámica de sistemas moleculares sometidos a cambios térmicos y volumétricos controlados. Entre los fenómenos estudiados destacan:

- Termalización del sistema desde un estado inicial fuera del equilibrio.
- Cristalización del argón mediante un proceso de enfriamiento lento y compresión.
- Evolución de propiedades como la temperatura, la presión, la energía potencial y la función de distribución radial \( g(r) \).

## Referencias fundamentales

- D. Frenkel, B. Smit, *Understanding Molecular Simulation*, Academic Press (2002).
- M. P. Allen, D. J. Tildesley, *Computer Simulation of Liquids*, Oxford University Press (1987).
- M. Tuckerman, *Statistical Mechanics: Theory and Molecular Simulation*, Oxford University Press (2010).

## Cómo ejecutar

Los códigos están escritos en Fortran 90/95. Se recomienda utilizar `gfortran`:

```bash
gfortran -ffast-math -O2 -o termostato termostato.f95
./termostato
