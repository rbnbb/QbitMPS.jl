# This Julia script creates documentation for project
# You should run this from the project root directory
# Ensure Documenter package is installed in the General registry
# run by typing $julia doc/make.jl
push!(LOAD_PATH, pwd())

using Documenter, QbitMPS

makedocs(; sitename = "QbitMPS")
