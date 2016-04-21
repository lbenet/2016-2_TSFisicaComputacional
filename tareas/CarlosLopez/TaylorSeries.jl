"""
Modulo de Taylor Series

Taylor(Array[])
Base: +, -, *, /, ^

Autor: Carlos López
22/03/2016
"""
__precompile__()

module TS

type Taylor{T<:Number} <: Number
    coffs :: Array{T,1}
    order :: Int

    function Taylor(coffs::Array{T,1}, order::Int)
        lencoef = length(coffs)
        order = max(order, lencoef-1)
        order == lencoef-1 && return new(coffs, order)
        resize!(coffs, order+1)
        for i = lencoef+1:order+1
            coffs[i] = zero(T)
        end
        new(coffs, order)
    end
end
Taylor{T<:Number}(x::Taylor{T}, order::Int) = Taylor{T}(x.coffs, order)
Taylor{T<:Number}(x::Taylor{T}) = x
Taylor{T<:Number}(coffs::Array{T,1}, order::Int) = Taylor{T}(coffs, order)
Taylor{T<:Number}(coffs::Array{T,1}) = Taylor{T}(coffs, length(coffs)-1)
Taylor{T<:Number}(x::T, order::Int) = Taylor{T}([x], order)
Taylor{T<:Number}(x::T) = Taylor{T}([x], 0)
Taylor{T<:Irrational}(coffs::Array{T,1}) =  Taylor(convert(Array{Float64,1}, coffs))
Taylor{T<:Irrational}(x::T) = Taylor(Array{Float64}([x]), 0)

import Base: print,show
function show(io::IO,x::Taylor)
    
    const simb = [c for c in "₀₁₂₃₄₅₆₇₈₉"]
    function sup(n::Int)
        dig = reverse(digits(n))
        join([simb[i+1] for i in dig])
    end
    
    n=length(x.coffs)-1
    
    for i = 1:n
        print(io,x.coffs[i],"ₜ",sup(i-1)," + ")
    end
    print(io,x.coffs[n+1],"ₜ",sup(n))
end
print(io::IO,x::Taylor)=show(io,x)

import Base: +, -, *, /, ^, == 
for f in (:+, :-)
    @eval begin
        function $f(A::Taylor,B::Taylor)
            a=zeros(Number,max(length(A.coffs),length(B.coffs)))
            for i=1:length(A.coffs)
                a[i]=A.coffs[i]
            end
            for i=1:length(B.coffs)
                a[i]=$f(a[i],B.coffs[i])
            end
            Taylor(a)
        end
        $f(x::Taylor,y::Number)=$f(x,Taylor(y))
        $f(y::Number,x::Taylor)=$f(x,Taylor(y))
    end
end
function *(A::Taylor,B::Taylor)
    a=zeros(Number,length(A.coffs)+length(B.coffs)-1)
    for k=1:length(a)
        for i=1:k
            if(i<=length(A.coffs) && k-i+1<=length(B.coffs))
                a[k]+=A.coffs[i]*B.coffs[k-i+1]
            end
        end
    end
    Taylor(a)
end
*(x::Taylor,y::Number)=x*Taylor(y)
*(y::Bool,x::Taylor)=x*Taylor(y)
*(y::Number,x::Taylor)=x*Taylor(y)

function degree(A::Taylor)
    l=length(A.coffs)
    for i=0:l-1
        if(A.coffs[l]!=0)
            break
        end
        l-=1
    end
    l-1
end

function shift(A::Taylor,d=1)
    if(d>0)
        a=zeros(Number,length(A.coffs)+d)
        for i=1:length(A.coffs)
            a[i+d]=A.coffs[i]
        end
        Taylor(a) 
    elseif(d<0 && length(A.coffs)+d>0)
        a=zeros(Number,length(A.coffs)+d)
        for i=1:length(a)
            a[i]=A.coffs[i-d]
        end
        Taylor(a)
    elseif(length(A.coffs)+d<1)
        Taylor(0)
    elseif(d==0)
        A
    end
end
    
function /(A::Taylor,B::Taylor)
    if(degree(B)==-1)
        error("division by zero")
    end
    if(degree(A)<degree(B))
        Taylor(0)
    else
        q=Taylor(zeros(Number,degree(A)-degree(B)+1))
        while(degree(A)>=degree(B))
            b=shift(B,degree(A)-degree(B))
            q.coffs[degree(A)-degree(B)+1]=A.coffs[degree(A)+1]/b.coffs[degree(b)+1]
            b=b*q.coffs[degree(A)-degree(B)+1]
            A=A-b
        end
        q
    end
end
/(x::Taylor,y::Number)=x/Taylor(y)
/(y::Number,x::Taylor)=x/Taylor(y)

function ==(A::Taylor,B::Taylor)
    R=true;
    if(length(A.coffs)!=length(B.coffs))
        a=zeros(Number,max(length(A.coffs),length(B.coffs)))
        b=zeros(Number,max(length(A.coffs),length(B.coffs)))
        for i=1:length(A.coffs)
            a[i]=A.coffs[i]
        end
        for i=1:length(B.coffs)
            b[i]=B.coffs[i]
        end
        for i=1:length(a)
            if(a[i]!=b[i])
                    R=false
            end
        end
    else
        for i=1:length(A.coffs)
            if(A.coffs[i]!=B.coffs[i])
                    R=false
            end
        end
    end
    R
end

import Base: exp,log,^,sin,cos
function exp(A::Taylor)
    a=zeros(Number,length(A.coffs))
    a[1]=exp(A.coffs[1])
    for k=2:length(a)
        for j=1:k-1
            a[k]+=(k-j)*A.coffs[k-j+1]*a[j]
        end
        a[k]=a[k]/(k-1)
    end
    Taylor(a)
end
function log(A::Taylor)
    d=zeros(Number,length(A.coffs)-1)
    for k=1:length(d)
        d[k]=A.coffs[k+1]*(k)
    end
    a=Taylor(d)/A
    r=zeros(Number,length(A.coffs))
    for i=2:length(a.coffs-1)
        r[i]=a.coffs[i-1]/(i-1)
    end
    r[1]=log(A.coffs[1])
    Taylor(r)
end
function ^(A::Taylor,x::Integer)
    if(x==0)
        Taylor(1)
    elseif(x>0)
        shift(A,x-1)
    else
        shift(A,x)
    end
end
^(A::Taylor,x::Rational)=exp(x*log(A))
^(A::Taylor,x::Real)=exp(x*log(A))

function sin(A::Taylor)
    a=zeros(Number,length(A.coffs))
    r=Taylor(a)
    for n=0:length(a)+1
        r+=(((-1)^n)/(factorial((2*n)+1)))*(A^((2*n)+1))
    end
    r
end
function cos(A::Taylor)
    a=zeros(Number,length(A.coffs))
    r=Taylor(a)
    for n=0:length(a)+1
        r+=(((-1)^n)/(factorial(2*n)))*(A^(2*n))
    end
    r
end

export Taylor

end