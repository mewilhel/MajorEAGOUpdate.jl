# Currently no refinement on min/max... will add later

"""
Reverse max
"""
max_rev!(a::MC, b::MC, c::MC) = ()
max_rev(a,b,c) = max_rev(promote(a,b,c)...)

"""
Reverse min
"""
min_rev!(a::MC, b::MC, c::MC) = ()
min_rev(a,b,c) = min_rev(promote(a,b,c)...)
