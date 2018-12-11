```
    EAGO.Fathom!

Selects and deletes nodes from stack with lower bounds greater than global
upper bound.
```
function Fathom!(d::Optimizer)
  if ~isempty(d.Stack)
    # Find the lowest upper bound without delecting
    MinKey, MinNode = first(d.Stack)
    for (Key, Node) in d.Stack
      if Node.UpperBound < d.GlobalUpperBound
        MinKey = Key
        MinNode = Node
        d.GlobalUpperBound = Node.UpperBound
      end
    end
    # Deletes all nodes with upper bound greater than minimum
    for (Key, Node) in d.Stack
      if Node.LowerBound > d.GlobalUpperBound
          delete!(d.Stack, Key)
      end
    end
  end
  d.CurrentNodeCount = length(d.Stack)
end
