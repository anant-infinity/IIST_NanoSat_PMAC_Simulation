

for k= 1:8016868
    k
    beta_value(k) = atan(norm(cross([B_body_x(k);B_body_y(k);B_body_z(k)],[B_eci_x(k);B_eci_y(k);B_eci_z(k)])/dot([B_body_x(k);B_body_y(k);B_body_z(k)],[B_eci_x(k);B_eci_y(k);B_eci_z(k)])));
end

plot(beta_value)