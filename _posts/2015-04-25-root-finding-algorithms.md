
## Root-Finding Algorithms

Here we discuss some of the numerical methods for finding the roots of the equations of the form **f(x) = 0**, where **'f'** is a function. We can find the roots in closed form for only some limited class of equations. Most of the time it is difficult to find roots of transedental equations. Sometimes, we can use graphical techniques to have rough idea of roots for such equations.

There are a variety of techniques to estimate the roots of the equation depending on the nature of the function involved.
Broadly we classify the different methods into two categories:
* **Root Bracketing Methods**
* **Iterative Methods**

Bracketing method is used when the root **'c'** of **f(x) = 0** is known to lie in the interval **(a,b)**. A possible guess of interval **(a,b)** for a given equation can be made by checking the signs of **f(a)** and **f(b)**. If signs are opposite then there is a possibility of root of the equation.

**Bisection Method** is one of the simplest numerical techniques based on root bracketing methods.
This method works if function **'f'** is continuous and **f(a)\*f(b) < 0**.
We find the sign of **f(x)** at the mid of the interval **(a,b)** i.e., **c = (a+b)/2** and check the sign of **f(c)**. If it matches with the sign of **f(a)** the sign of then the interval **(a,b)** is replaced with **(c,b)** else with **(a,b)**. This process is repeated and it can be stopped by one of the following mechanisms:
* the absolute value of the difference between two consecutive estimates of mid-points of the interval satisfies some tolerance limit.
* the absolute value of the function **'f'** at the mid-point of the interval satisfies some tolerance limit.

Here is a sample python code to estimate the one of the roots of **$$ 2x^4-7x^3-15x^2+34x+40 = 0 $$** using Bisection method. Clearly, **$ f(x) = 2x^4-7x^3-15x^2+34x+40 $** and **(a,b)** can be taken as (0,3).


```python
from __future__ import  division
import numpy as np

def f(x):
	return 2*x**4-7*x**3-15*x**2+34*x+40

a = input ("Enter the lower bound of the interval in which root may lie:")
b = input ("Enter the upper bound of the interval in which root may lie:")

c = (a+b)/2

if(f(a)*f(b) > 0) :
	print a, "  and ",b," are not appropriate for bisection method."
else: 
	while ((b-a) > 0.001):
		if(f(a)*f(c) < 0):
			b=c
			c=(a+b)/2
		else :
			a=c
			c=(a+b)/2
	print c," is the approximate root."
```

    Enter the lower bound of the interval in which root may lie:0
    Enter the upper bound of the interval in which root may lie:3
    2.50012207031  is the approximate root.


Let us plot the above function  **$ f(x) = 2x^4-7x^3-15x^2+34x+40 $** to see the location of roots.


```python
%matplotlib inline
import matplotlib.pyplot as plt

x=np.linspace(-3,5,100)
plt.xlabel('x')
plt.ylabel(r'$f(x)$')
plt.plot(x,2*x**4-7*x**3-15*x**2+34*x+40)
plt.plot(x,x*0)
plt.axis([-3,5,-30,60])
plt.show()
```


![png](Untitled_files/Untitled_7_0.png)


There is an another technique known as **Regula-Falsi Method** under Root Bracketing Method which has better convergence than Bisection Method. Here to choose the point **'c'**, we also consider the role of function **f(x)** along with chosen interval **(a,b)**,

$$c = b-\frac{(b-a)*f(b)}{f(b)-f(a)} = \frac{a*f(b)-b*f(a)}{f(b)-f(a)}$$

Now let us discuss some Iterative Methods for numerical approximation of roots. Here, we start with a guess of the root **'$x_0$'** and using some recursion relation we get better and better approximations of the root. The method works if the approximate roots generated recursively converge.

**Newton's Method** or **Newton-Raphson Method** is one of such methods. Here to find the next approximation of root, stating with **$x_0$**, we calculate the x-intercept of tangent line to the curve **f(x)** at the point corresponding to **$x_0$**.
$$x_1 = x_0 - \frac{f(x_0)}{f'(x_0)}$$
In general the recursion relation can be written as,
$$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$
This recursion can be stopped by
* fixing apriori the total number of iterations
* setting a tolerance limit for $|x_{n+1}-x_n|$


This method works fine if **f(x)** is differentiable and convergence is better if the initial guess $x_0$ is close enough to actual unknown root.

Here is a sample python code to estimate the root of **$x- exp(-x^2) = 0 $**, using Newton's method. Here, **$ f(x) = x- exp(-x^2) $**, **$ f'(x) = 1+2x* exp(-x^2) $** and the guess $x_0$ can be taken as 1.0


```python
from __future__ import  division
import numpy as np

def f(x):
	return x-np.exp(-x**2)

def f_(x):
	return 1+2*x*np.exp(-x**2)     # or one can estimate derivative with 
                                   #     (f(x+h)-f(x-h))/(2*h)
 
x1 = input ("Enter the guess value of root:")

x2 = x1-(f(x1)/f_(x1))

while (abs(x2-x1)> 0.001):
	x1=x2
	x2 = x1-(f(x1)/f_(x1))
	
print x2, "is the approximate root."
```

    Enter the guess value of root:1.0
    0.652918640437 is the approximate root.


The plot of $f(x)=x-exp(-x^2)$ looks like


```python
%matplotlib nbagg
import matplotlib.pyplot as plt
import numpy as np

x=np.linspace(-2,5,1000)
plt.xlabel('x')
plt.ylabel(r'$f(x)$')
plt.plot(x,x-np.exp(-x**2))
plt.plot(x,x*0)
plt.show()
```


    <IPython.core.display.Javascript object>



<img src="data:image/png;base64,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" width="432">


There are many other iterative methods of root approximation like **Secant Method, Muller's Method**, which may be more efficient than Newton's Method. In Secant method, we use method of finite difference to evaluate the derivative (which can be helpful if $f(x)$ has complicated form).
