# [Particle Filter](https://www.wikiwand.com/en/Particle_filter)

## Introduction

[A good video explaining Particle Filter clearly and vividly](https://www.youtube.com/watch?v=NrzmH_yerBU)

## Demos

|demo1| demo2|
|:---:|:---:|
| Sampling globally/locally | Heuristic policy|
|||

## Remarks

You may be confused why sometimes the sum of particle probabilities is zero.

It can happen in two main situations:

1. When lidar scan across a concave corner, the true measurement suddenly changes.

Before the change, the particles which already pass the corner vanish. And after the change, all particles that still not pass the corner are "not trusted".

2. When lidar starts to scan the uncertain bound, or end scanning.

Before scaning the uncertain bound, the estimation of the robot position always converges. So do the estimation of the wall. Consequently, when lidar starts to scan the uncertain bound, the true measurement differs a lot from the estimated one.

Solution: 

I offer 4 remedy policies to tackle this problem. You can check them in the code.

## File Structure

| file name| notes|
|:---:|:---:|
| Estimator.m | main implementation |
| run.m | simulation & visulization |
| ...Consts | config params here |

### How to play? 

run "run.m" in MATLAB






