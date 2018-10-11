```
    EAGO.FindLowerBound

Selects node with the lowest lower bound in d.Stack.
```
function FindLowerBound(d::Optimizer)
  stack = d.Stack
  minkey, minnode = next(stack, start(stack))[1]
  minvalue = minnode.LowerBound
  for (key, node) in stack
    if node.LowerBound < minvalue
      minkey = key
      minnode = node
      minvalue = node.LowerBound
    end
  end
  minvalue
end
