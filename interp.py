import numpy as np
import matplotlib.pyplot as plt

#y_data = np.array([2.68e6,3.88e6,4.24e6]) #RUN: 62,63,64
y_data = np.array([2.66e6,3.89e6,4.21e6]) #RUN: 67,68,69

x_data = np.array([4.5,6.5,7])

linear_model=np.polyfit(x_data,y_data,1)
linear_model_fn=np.poly1d(linear_model)
z = 3.5
f_z = linear_model_fn(z)

print("Interpolation value for %.2fV (OV) = %.2f"%(z,f_z))

x = np.linspace(2,8)
plt.plot(x, linear_model_fn(x), label = linear_model_fn, color = "k")
plt.scatter(x_data, y_data, label = "data", color = "b")
plt.scatter(z,f_z,label = "Interpolated value f(%.2f)=%.2f"%(z,f_z), color = "r")
plt.plot(x,f_z*np.ones(np.size(x)),color = "r", ls = "--")

plt.title("Gain interpolation SC (RUN 67-69)")
plt.xlabel("OV in [V]");plt.ylabel("Gain in [e-]")
plt.legend()
plt.show()  