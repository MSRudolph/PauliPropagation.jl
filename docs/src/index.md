# API


## Index

```@index
```


## Functions

```@docs
propagate(circ, pstr::PauliString{TT,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:PauliStringType}
```

```@docs
propagate(circ, psum::PauliSum{TT,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:PauliStringType}
```

```@docs
propagate!(circ, psum::PauliSum{TT,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:PauliStringType}
```
