export set_config, get_config

__config=Dict{String,Any}() 

function get_config(key::String)
	return __config[key]
end

function get_config(key::String, default_value)
	return get(__config, key, default_value)
end

function set_config(key::String,value)
	__config[key]=value
end

using DataStructures
global calldict = DefaultDict{String}{Int}(0) 
