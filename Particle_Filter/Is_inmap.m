function result = Is_inmap(x,y)
    cond1 = x<=0.5 && y<=0.5;
    cond2 = x<=1 && y>=2.5;
    cond3 = x>=1 && x<=1.25 && y>=-x+3.5;
    cond4 = x>=1.25 && x<=2 && y>=x+1;
    cond5 = x>=2 && y>= -x+5;
    cond6 = x>=2.5 && y<= x-1;
    if cond1||cond2||cond3||cond4||cond5||cond6
        result = false;
    else
        result = true;
    end
end