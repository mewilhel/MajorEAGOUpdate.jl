```
    EAGO.FindLowerBound

Selects node with the lowest lower bound in d.Stack.
```
function FindLowerBound(d::Optimizer)
  minkey, minnode = first(d.Stack)
  minvalue = minnode.LowerBound
  for (key, node) in d.Stack
    if node.LowerBound < minvalue
      minkey = key
      minnode = node
      minvalue = node.LowerBound
    end
  end
  minvalue
end
