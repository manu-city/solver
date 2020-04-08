# Time Marching

## Dependencies
- mesh
- solution
- CFL
- sigma
- nonDimParams

## Initialisation of variables

    u = solution.u;
    v = solution.u;

    dx = mesh.dx;
    dy = mesh.dy;

    nu = nonDimParams.v_;
    
## dt calculation

    min_1 = min(dx./abs(u), [], 'all');
    min_2 = min(dy./abs(v), [], 'all');
    min_3 = min(dx^2/nu);

    dt = min(CFL * min(min_1, min_2), sigma * min_3);
