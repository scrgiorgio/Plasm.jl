export set_config, get_config

__config=Dict()


function get_config(key)
  return __config[key]
end

function get_config(key, default_value)
  return get(__config, key, default_value)
end

function set_config(key,value)
  __config[key]=value
end

using DataStructures
global calldict = DefaultDict{String}{Int}(0) 
