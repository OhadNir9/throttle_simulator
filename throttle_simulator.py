import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import numpy as np
import tkinter as tk
from tkinter import ttk

"""
radio_min = 1000
radio_max = 2000
deadzone = 30
thr_dz = 100
high_in = 1000
pilot_speed_dn = 12
pilot_speed_up = 8


radio_trim_low = radio_min + deadzone
"""
HIGH_IN = 1000


def constrain_int(input, min, max):
    if input < min:
        return min
    if input > max:
        return max
    return input

def get_throttle_mid(radio_min, radio_max, radio_trim_low):
    return HIGH_IN * (((radio_min + radio_max)/2) - radio_trim_low) / (radio_max - radio_trim_low)


def calculate_pwm_speed_series(radio_min, radio_max, deadzone, thr_dz, pilot_speed_dn, pilot_speed_up):
    pwm_speed_vals = []
    pwms = []
    speeds = []
    radio_trim_low = radio_min + deadzone
    for pwm_input in range (1000, 2100, 10):
        ## pwm_to_range_dz ##
        r_in = constrain_int(pwm_input, radio_min, radio_max)


        if r_in > radio_trim_low:
            control_in = (((HIGH_IN) * (r_in - radio_trim_low)) / (radio_max - radio_trim_low))
        else:
            control_in = 0

        ## get_pilot_desired_climb_rate ##
        throttle_control = constrain_int(control_in, 0, 1000)
        mid_stick = get_throttle_mid(radio_min, radio_max, radio_trim_low)
        deadband_top = mid_stick + thr_dz
        deadband_bottom = mid_stick - thr_dz
        # check throttle is above, below or in the deadband
        if (throttle_control < deadband_bottom):
            # below the deadband
            desired_rate = pilot_speed_dn * (throttle_control-deadband_bottom) / deadband_bottom
        elif (throttle_control > deadband_top):
            # above the deadband
            desired_rate = pilot_speed_up * (throttle_control-deadband_top) / (1000-deadband_top)
        else:
            # must be in the deadband
            desired_rate = 0

        pwm_speed_vals.append((pwm_input, desired_rate))
        pwms.append(pwm_input)
        speeds.append(desired_rate)
    #return pwm_speed_vals
    return pwms, speeds

initial_values = {'radio_min': 1000, 'radio_max': 2000, 'deadzone': 30, 'thr_dz': 100, 'pilot_speed_dn': 12, 'pilot_speed_up': 8}
# Create initial plot
fig, ax = plt.subplots(figsize=(15,15))
plt.subplots_adjust(bottom=0.4)
x, y = calculate_pwm_speed_series(**initial_values)
ax.set_ylim([-15, 15])
ax.set_xlim([1000, 2000])
line, = ax.plot(x, y)

ax.set_title("Throttle Simulator!", {'fontsize': 20})

# Create sliders for each constant
#slider_ax = plt.axes([0.1, 0.01, 0.65, 0.03])
#slider_a = Slider(slider_ax, 'a', -10, 10, valinit=initial_values['a'])
slider_radio_min = Slider(plt.axes([0.17, 0.01, 0.65, 0.03]), 'radio_min', 1000, 1150, valinit=initial_values['radio_min'])
slider_radio_max = Slider(plt.axes([0.17, 0.06, 0.65, 0.03]), 'radio_max', 1850, 2000, valinit=initial_values['radio_max'])
slider_deadzone = Slider(plt.axes([0.17, 0.11, 0.65, 0.03]), 'deadzone', 0, 400, valinit=initial_values['deadzone'])
slider_thr_dz = Slider(plt.axes([0.17, 0.16, 0.65, 0.03]), 'thr_dz', 0, 400, valinit=initial_values['thr_dz'])
slider_pilot_speed_dn = Slider(plt.axes([0.17, 0.21, 0.65, 0.03]), 'pilot_speed_dn', 2, 15, valinit=initial_values['pilot_speed_dn'])
slider_pilot_speed_up = Slider(plt.axes([0.17, 0.26, 0.65, 0.03]), 'pilot_speed_up', 2, 15, valinit=initial_values['pilot_speed_up'])

joystick_pwm_limits = [None, None] # Initial values
def calc_joystick_pwm_limits(label):
    global joystick_pwm_limits
    if label == "Slow":
        joystick_pwm_limits = [1350, 1650]
    if label == "Medium":
        joystick_pwm_limits = [1200, 1800]
    if label == "Fast":
        joystick_pwm_limits = [1000, 2000]
    print(joystick_pwm_limits)
    
radio = RadioButtons(plt.axes([0.9, 0.9, 0.1, 0.1]), ('Slow', 'Medium', 'Fast'),
                     active=0)
radio.on_clicked(calc_joystick_pwm_limits)
calc_joystick_pwm_limits(radio.value_selected)

intersections_y = [None, None]
intersect_lines_y = [None, None]
intersect_lines_x = [None, None]
intersect_annotations = [None, None]

for idx, pwm in enumerate(joystick_pwm_limits):
    # Find intersection point of vertical line and the graph
    intersections_y[idx] = np.interp(pwm, x, y)
    # Add horizontal line at the intersection point
    intersect_lines_y[idx] = ax.axhline(y=intersections_y[idx], color='g', linestyle='--')
    # Add vertical line at x=
    intersect_lines_x[idx] = ax.axvline(x=pwm, color='r', linestyle='--', ymin=0, ymax=((intersections_y[idx]+15)/30))

    if idx == 0:
        lim_name = "min"
    else:
        lim_name = "max"
    intersect_annotations[idx] = ax.annotate(f'V{lim_name}[m/s]={intersections_y[idx]:.2f}', xy=(1000, intersections_y[idx]), xytext=(5, 0),
                textcoords='offset points', color='g', fontsize=16)

# Function to update the plot when sliders change
def update(val):
    global intersect_lines_y
    global intersect_lines_x
    global intersect_annotations
    global joystick_pwm_limits

    a = slider_radio_min.val
    b = slider_radio_max.val
    c = slider_deadzone.val
    d = slider_thr_dz.val
    e = slider_pilot_speed_dn.val
    f = slider_pilot_speed_up.val
    x, y = calculate_pwm_speed_series(a, b, c, d, e, f)
    line.set_xdata(x)
    line.set_ydata(y)

    for idx, pwm in enumerate(joystick_pwm_limits):
        intersect_lines_y[idx].remove()
        intersect_lines_x[idx].remove()
        # Find intersection point of vertical line and the graph
        intersections_y[idx] = np.interp(pwm, x, y)
        # Add horizontal line at the intersection point
        intersect_lines_y[idx] = ax.axhline(y=intersections_y[idx], color='g', linestyle='--')
        # Add vertical line at x
        intersect_lines_x[idx] = ax.axvline(x=pwm, color='r', linestyle='--', ymin=0, ymax=((intersections_y[idx]+15)/30))
    
        # Annotate the horizontal line with its y-value on the y-axis
        intersect_annotations[idx].remove()
        if idx == 0:
            lim_name = "min"
        else:
            lim_name = "max"
        intersect_annotations[idx] = ax.annotate(f'V{lim_name}[m/s]={intersections_y[idx]:.2f}', xy=(1000, intersections_y[idx]), xytext=(5, 0),
                textcoords='offset points', color='g', fontsize=20)

    fig.canvas.draw_idle()

# Connect sliders to update function
slider_radio_min.on_changed(update)
slider_radio_max.on_changed(update)
slider_deadzone.on_changed(update)
slider_thr_dz.on_changed(update)
slider_pilot_speed_dn.on_changed(update)
slider_pilot_speed_up.on_changed(update)
radio.on_clicked(update)

plt.show()