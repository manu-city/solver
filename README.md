# Boundary of Fluid method

## Dependencies
- mesh
- nonDimParams
- dimParams
- N
- M
- prediction

## Working
Grid points in x between ribs
`
spacing = clear/dx;
`

Grid points in y rib occupies
`o_height = o_y/dy;`

Domain Matrix
`VOF = ones(N,M);`

Allocation of 0's in domain matrix

  
      for i = N-(o_height - 1):N
        for n = 1:o_num
            switch n
                case 1 
                    for j = spacing/2 + 1:(spacing/2 + 1) + o_length
                        VOF(i,j) = 0;    
                    end
                otherwise 
                    for j = spacing * (n) + (2*n-3):spacing * (n) + (2*n-3) + o_length
                        VOF(i,j) = 0;
                    end
            end 
        end
      end

Updating prediction with null points

    prediction.predU = u .* VOF;
    prediction.predV = v .* VOF;
