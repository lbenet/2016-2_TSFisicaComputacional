include("AutomDiff1.jl")
using Base.test
using AD

A=(0,1,7,0.0,1.0,7.0,pi,-1,-7)
for (a in A)
    for(b in A)
        d=Dual(a,b);
        c=xdual(a)
        if(a!=0)
            @test (d+d-d*d/d)^1==d
            @test (d+a-a*a/a)^1==d
            @test (c+d-d*d/d)^1==c
        else
            @test ((d+d-d)*d)^1==d*d
            @test ((d+a-a)*a)^1==d*a
            @test ((c+d-d)*d)^1==c*d
        end
    end
end

derivadas = (
    [:exp,:(exp(x.fun))],
    [:sqrt,:(1/(2*sqrt(x.fun)))],
    [:log,:(1/x.fun)],
    [:sin,:(cos(x.fun))],
    [:cos,:(-sin(x.fun))],
    [:tan,:(sec(x.fun)^2)],
    [:cot,:(-csc(x.fun)^2)],
    [:sec,:(sec(x.fun)*tan(x.fun))],
    [:csc,:(-csc(x.fun)*cot(x.fun))],
    [:sinh,:(cosh(x.fun))],
    [:cosh,:(sinh(x.fun))],
    [:tanh,:(sec(x.fun)^2)],
    [:coth,:(-csch(x.fun)^2)],
    [:sech,:(sech(x.fun)*tanh(x.fun))],
    [:csch,:(-csch(x.fun)*coth(x.fun))],
    [:asin,:(1/sqrt(1-(x.fun^2)))],
    [:acos,:(-1/sqrt(1-(x.fun^2)))],
    [:atan,:(1/sqrt(1+(x.fun^2)))],
    [:acot,:(-1/(1+(x.fun^2)))],
    [:asec,:(1/(abs(x.fun)*sqrt((x.fun^2)-1)))],
    [:acsc,:(-1/(abs(x.fun)*sqrt((x.fun^2)-1)))],
    [:asinh,:(1/sqrt((x.fun^2)+1))],
    [:acosh,:(1/sqrt((x.fun^2)-1))],
    [:atanh,:(1/sqrt(1-(x.fun^2)))],
    [:acoth,:(1/sqrt(1-(x.fun^2)))],
    [:asech,:(1/(1-(x.fun^2)))],
    [:acsch,:(-1/(abs(x.fun)*sqrt((x.fun^2)+1)))],
    [:sinc,:(cosc(x.fun))],
)

using Base.Test
A=(1,7,0.0,1.0,7.0,pi,-1,-7)
x=0
for (a in A)
    for(b in A)
        x=Dual(a,b)
        for op in derivadas
            da=false;
            db=false;
            try
                da=eval(quote$(op[1])end)(x)
            end
            try
                db=Dual(eval(quote$(op[1])end)(a) , b*(eval(quote$(op[2])end)))
            end
            print("da ",da,"\ndb ",db,"\n")
            if(da!=false && db!=false && !isnan(db.fun) && !isnan(db.der))
                print ("tested: ",op[1],"(",x,") = ",da,"\n",)
                @test da==db
            end
        end
    end
end