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
```
The most widly used function are "and"/"or"(intersection/union) operation.
```python
"|" means or, "&" means and.
(P.open(1, 11) | P.closed(0, 1) | P.closed(20, 21)) # [0, 11) | [20, 21]
P.open(1, 11) & P.closed(0, 1) # () minus P.empty()
P.closed(1, 11) & P.closed(0, 1) # [1], note that I use P.closed for (1, 11)
```
