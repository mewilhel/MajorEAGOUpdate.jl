```
    EAGO.Fathom!

Selects and deletes nodes from stack with lower bounds greater than global
upper bound.
```
function Fathom!(d::Optimizer)
  # Find the lowest upper bound without delecting
  Stack = d.Stack
  MinKey, MinNode = next(Stack, start(Stack))[1]
  GlobalUpperBound = MinNode.UpperBound
  for (Key, Node) in Stack
    if Node.UpperBound < MinValue
      MinKey = Key
      MinNode = Node
      GlobalUpperBound = Node.UpperBound
    end
  end
  # Deletes all nodes with upper bound greater than minimum
  for (Key, Node) in Stack
    if node.LowerBound > GlobalUpperBound
        delete!(Stack)
    end
  end
end
