% DEPENDENT ON THE MAP CONFIGURATION!!
% HARD CODING HERE!!

% how to know whether a ray and a line segment intersect?
% the algorithm is adapted from the website below
% https://stackoverflow.com/questions/14307158/
function distance = get_distance(xr,yr,phi,kappa,contour)
    contour(8,1) = kappa;
    contour(9,1) = kappa;
    contour = [contour; contour(1,:)];
    distance = 99999;
    for seg = 1:10
        % (x1,y1) --(x,y)-- (x2,y2)
        x1 = contour(seg,1); y1 = contour(seg,2);
        x2 = contour(seg+1,1); y2 = contour(seg+1,2);
        % (x,y) = (x1,y1) + u(x2-x1,y2-y1)= (xr,yr) + t(cos(phi),sin(phi))
        A = [x2 - x1, -cos(phi);
             y2 - y1, -sin(phi)];
        b = [xr - x1;
             yr - y1];
        if rank(A)~=rank([A,b]) % no solution
        elseif rank(A) < 2 && dot([x1-xr; y1-yr],[cos(phi);sin(phi)])>0% inf solutions => collinear
            distance = min([distance,norm([x1-xr; y1-yr]),norm([x2-xr; y2-yr])]);
        elseif rank(A) == 2 % unique solution
            C = A\b;
            if C(1)>=0 && C(1)<=1 && C(2)>=0
                xc = x1 + C(1)*(x2 - x1);
                yc = y1 + C(1)*(y2 - y1);
                distance = min(distance,norm([xc-xr; yc-yr]));
            end
        end
    end
end