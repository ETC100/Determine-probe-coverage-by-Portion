# Determine-probe-coverage-by-Portion
# What is portion?
Portion is a python package for set operation. If you use it well, you can even achieve a "pseudo bedtools" in your python program.  
The first thing we need to do is to know how to use portion, and then, we will learn how to determine the probe region for disease monitoring.  
Let's install the portion now.
```
pip install portion
```
and import it
```python
import portion as P # import portion
interval = P.closedopen(2, 4) #(CLOSED, 2, 4, OPEN), four args here, two for open/closed, two for number recording.
```
The document for portion is here. I will help you to remember the functions in portion package.  
https://pypi.org/project/portion/  
```python
## "closed" on the left or right means closed interval on the left or right, as for "open", it's the same.
P.open(1, 2) # (1,2)
P.closed(1, 2) # [1,2]
P.openclosed(1, 2) # (1,2]
P.closedopen(1, 2) # {1,2)
P.singleton(1) # only to contain one number
P.empty() # create empty P object
interval = P.open(1, 2)
interval.empty # judge if the interval contain nothing, here return False P.open(0, 0) True
```
The most widly used function are "and"/"or"(intersection/union) operation.
```python
# "|" means or, "&" means and.
(P.open(1, 11) | P.closed(0, 1) | P.closed(20, 21)) # [0, 11) | [20, 21]
P.open(1, 11) & P.closed(0, 1) # () minus P.empty()
P.closed(1, 11) & P.closed(0, 1) # [1], note that I use P.closed for (1, 11)
~P.closed(0, 1) # reverse operation, (-inf,0) | (1,+inf) inf, P.inf
P.closed(1,2) in P.closed(0, 3) # True, the left interval surely in the right
```
A portion object may contain one or several intervals, we can use atomic to judge if there is only one in it.  
This is usefu
```python
interval = (P.open(1, 11) | P.closed(0, 1) | P.closed(20, 21))
interval.atomic # return False, here are two intervals
P.open(0, 1).atomic # return True
```
Acquire the smallest interval covering provided portion object
```python
(P.closed(0, 1) | P.open(2, 3)).enclosure #[0,3)
interval = (P.open(1, 11) | P.closed(0, 1) | P.closed(20, 21))
interval.enclosure # [0 ,21]
```
Modfiy portion object needs replace function, no exterior setter.
```python
i = P.closed(0, 2)
i.replace(P.OPEN, -1, 3, P.CLOSED) # (-1,3]
```
Portion can also work as the dictionary.
```python
## dictionary
d = P.IntervalDict()
d[P.closed(0, 3)] = 'banana'
d[4] = 'apple'
{[0,3]: 'banana', [4]: 'apple'}
d.find('apple') # 4

list(d.keys()) # [[0,2), [2,4]]
list(d.values()) # ['banana', 'orange']
list(d.items()) # [([0,2), 'banana'), ([2,4], 'orange')]
```
Output and Input
```python
s = P.to_string(P.closedopen(0, 1)) # '[0,1)', we can store this string in somewhere
interval = P.from_string(s) # read the string
interval # [0,1), back to portion object
```

Okay, let's start to design the probe for disease monitoring. Here we use pan-cancer as the exmaple.  
The first thing is to find out the regions we should monitor.  
Data we need:  
1. Known conductive regions, from reference or consensus
2. Your in-house regions, from your experience
