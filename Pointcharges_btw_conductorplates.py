import sympy as sp


# Define the variables and constants
z, r   = sp.symbols('z r')
mu, z0 = sp.symbols('mu z0', real = True, constant = True)

# Define the Green function explicitly
#f = sp.cos(z * r)
f = sp.sqrt( 1.0 + sp.tan(z*mu)**2 * sp.tanh(r*mu)**2 ) * ( \
        1.0/sp.sqrt( ( sp.tan(z*mu) - sp.tan(z0*mu) )**2 + sp.tanh(r * mu)**2 * \
                 ( 1.0 + sp.tan(z*mu) * sp.tan(z0*mu) )**2 ) - \
        1.0/sp.sqrt( ( 1.0 - sp.tan(z*mu) * sp.tan(z0*mu) )**2 + sp.tanh(r * mu)**2 * \
                 ( sp.tan(z*mu) + sp.tan(z0*mu) )**2 ) \
        )


# Compute the partial derivatives
df_dz = sp.diff(f, z)           # This will give us the field at the plates
df_dr = sp.diff(f, r)

# Define a specific value for z to evaluate the derivatives
z_val1 = sp.N(sp.pi)/(4.0*mu)   # Change this to the desired value
z_val2 = - sp.N(sp.pi)/(4.0*mu)

# Evaluate f at z=d/2 and z=-d/2
f_single1_arg = f.subs(z, z_val1)
f_single2_arg = f.subs(z, z_val2)

# Substitute z_values into the partial derivatives to get functions of r only
df_dz_single1_arg = df_dz.subs(z, z_val1)         # Gives the field at the plate in z=0 as a function of r
df_dr_single1_arg = df_dr.subs(z, z_val1)

df_dz_single2_arg = df_dz.subs(z, z_val2)         # Gives the field at the plate in z=D as a function of r
df_dr_single2_arg = df_dr.subs(z, z_val2)

# Print the results
print('')
print("Green function f(z, r, z0, 0):", f)
print('')
print("Partial derivative with respect to z:", sp.simplify(df_dz) )
print('')
print("Partial derivative with respect to r:", df_dr )
print('')
print('')
print(f"Partial derivative with respect to z evaluated at z = {z_val1}:", sp.simplify(df_dz_single1_arg) )
print('')
#print(f"Partial derivative with respect to r evaluated at z = {z_val1}:", df_dr_single1_arg)
#print('')
print('')
print(f"Partial derivative with respect to z evaluated at z = {z_val2}:", sp.simplify(df_dz_single2_arg) )
print('')
#print(f"Partial derivative with respect to r evaluated at z = {z_val2}:", df_dr_single2_arg)
#print('')


## Integrate the evaluated partial derivatives with respect to x
#integrated_df_dx = sp.integrate(df_dx_single_arg, x)
    
################################################################################################
################################################################################################
################################################################################################

## To make it a function of a single argument r
#f_x = sp.lambdify(x, df_dx_single1_arg, 'numpy')
#g_x = sp.lambdify(x, df_dy_single1_arg, 'numpy')
#
## Example of using the functions
#r_value = 1  # Change this to the desired value
#print(f"f_x({x_value}) = {f_x(x_value)}")
#print(f"g_x({x_value}) = {g_x(x_value)}")




#1.0*mu*(2.0*((1.0*(tan(mu*z0) - 1)*tan(mu*z0) + (tan(mu*z0) + 1.0)*tanh(mu*r)**2)*((tan(mu*z0) - 1.0)**2 + 1.0*(tan(mu*z0) + 1)**2*tanh(mu*r)**2)**(3/2) - (1.0*(tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1.0)**2*tanh(mu*r)**2)**(3/2)*(1.0*(tan(mu*z0) + 1)*tan(mu*z0)*tanh(mu*r)**2 - tan(mu*z0) + 1.0))*(tanh(mu*r)**2 + 1) + 2.0*(sqrt(1.0*(tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1.0)**2*tanh(mu*r)**2) - sqrt((tan(mu*z0) - 1.0)**2 + 1.0*(tan(mu*z0) + 1)**2*tanh(mu*r)**2))*(1.0*(tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1.0)**2*tanh(mu*r)**2)*((tan(mu*z0) - 1.0)**2 + 1.0*(tan(mu*z0) + 1)**2*tanh(mu*r)**2)*tanh(mu*r)**2)/((1.0*(tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1.0)**2*tanh(mu*r)**2)**(3/2)*((tan(mu*z0) - 1.0)**2 + 1.0*(tan(mu*z0) + 1)**2*tanh(mu*r)**2)**(3/2)*sqrt(tanh(mu*r)**2 + 1))
#=
#mu * (
#      2.0 * (
#             ( (tan(mu*z0) - 1) * tan(mu*z0) + (tan(mu*z0) + 1.0) * tanh(mu*r)**2 ) *
#             ( (tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1)**2 * tanh(mu*r)**2 )**(3/2)
#             -
#             ( (tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1)**2 * tanh(mu*r)**2 )**(3/2) *
#             ( (tan(mu*z0) + 1) * tan(mu*z0) * tanh(mu*r)**2 - tan(mu*z0) + 1.0 )
#            ) * (tanh(mu*r)**2 + 1)
#    + 2.0 * (
#             sqrt( (tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1)**2 * tanh(mu*r)**2 )
#           - sqrt( (tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1)**2 * tanh(mu*r)**2 )
#            )  * ( (tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1)**2 * tanh(mu*r)**2 ) *
#            ( (tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1)**2 * tanh(mu*r)**2 ) * tanh(mu*r)**2
#     )/(
#     ( (tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1)**2 * tanh(mu*r)**2 )**(3/2) *
#     ( (tan(mu*z0) - 1)**2 + (tan(mu*z0) + 1)**2 * tanh(mu*r)**2 )**(3/2) *
#       sqrt(tanh(mu*r)**2 + 1)
#       )

##################################################

