#ifndef SOIL_RHO_H
#define SOIL_RHO_H

float default_soil_vector_1[ ] = {0.081, 0.089, 0.098, 0.107, 0.116, 0.125, 0.133, 0.142, 0.151, 0.159, 0.168, 0.176, 0.185, 0.194, 0.202, 0.211, 0.219, 0.228, 0.236, 0.245, 0.254, 0.262, 0.271, 0.280, 0.289, 0.298, 0.308, 0.318, 0.328, 0.338, 0.349, 0.361, 0.373, 0.385, 0.398, 0.410, 0.421, 0.431, 0.441, 0.450, 0.459, 0.466, 0.474, 0.482, 0.489, 0.497, 0.504, 0.512, 0.519, 0.527, 0.534, 0.542, 0.549, 0.557, 0.565, 0.573, 0.581, 0.589, 0.598, 0.606, 0.615, 0.624, 0.633, 0.642, 0.651, 0.661, 0.670, 0.679, 0.688, 0.697, 0.706, 0.714, 0.723, 0.731, 0.739, 0.747, 0.754, 0.762, 0.769, 0.776, 0.782, 0.789, 0.794, 0.800, 0.805, 0.810, 0.814, 0.819, 0.822, 0.826, 0.829, 0.833, 0.836, 0.840, 0.844, 0.847, 0.851, 0.855, 0.859, 0.863, 0.867, 0.871, 0.876, 0.880, 0.884, 0.889, 0.894, 0.898, 0.903, 0.908, 0.913, 0.918, 0.924, 0.929, 0.935, 0.940, 0.946, 0.952, 0.958, 0.963, 0.969, 0.975, 0.980, 0.986, 0.992, 0.997, 1.003, 1.008, 1.014, 1.019, 1.025, 1.030, 1.035, 1.041, 1.046, 1.051, 1.056, 1.061, 1.065, 1.069, 1.074, 1.078, 1.081, 1.085, 1.088, 1.092, 1.095, 1.098, 1.100, 1.103, 1.106, 1.108, 1.111, 1.114, 1.116, 1.119, 1.122, 1.125, 1.128, 1.131, 1.133, 1.136, 1.139, 1.142, 1.146, 1.149, 1.153, 1.156, 1.160, 1.164, 1.168, 1.173, 1.177, 1.181, 1.184, 1.187, 1.189, 1.191, 1.191, 1.191, 1.191, 1.190, 1.188, 1.186, 1.183, 1.179, 1.174, 1.170, 1.164, 1.158, 1.151, 1.143, 1.136, 1.127, 1.117, 1.101, 1.081, 1.055, 1.020, 0.979, 0.935, 0.889, 0.845, 0.811, 0.788, 0.773, 0.762, 0.754, 0.748, 0.744, 0.743, 0.743, 0.746, 0.750, 0.757, 0.766, 0.777, 0.789, 0.803, 0.817, 0.832, 0.848, 0.864, 0.881, 0.898, 0.914, 0.930, 0.945, 0.960, 0.974, 0.988, 1.001, 1.014, 1.026, 1.037, 1.048, 1.058, 1.068, 1.077, 1.085, 1.093, 1.100, 1.107, 1.113, 1.119, 1.125, 1.130, 1.135, 1.139, 1.143, 1.147, 1.150, 1.153, 1.155, 1.157, 1.159, 1.161, 1.162, 1.162, 1.163, 1.163, 1.162, 1.161, 1.159, 1.157, 1.154, 1.151, 1.148, 1.144, 1.139, 1.134, 1.129, 1.122, 1.116, 1.109, 1.101, 1.094, 1.088, 1.082, 1.078, 1.074, 1.071, 1.066, 1.061, 1.056, 1.050, 1.044, 1.038, 1.030, 1.019, 1.004, 0.985, 0.962, 0.934, 0.895, 0.843, 0.780, 0.707, 0.625, 0.551, 0.495, 0.456, 0.422, 0.389, 0.360, 0.338, 0.326, 0.323, 0.323, 0.325, 0.329, 0.334, 0.340, 0.347, 0.355, 0.364, 0.374, 0.385, 0.397, 0.408, 0.420, 0.433, 0.446, 0.459, 0.472, 0.486, 0.500, 0.514, 0.528, 0.542, 0.555, 0.569, 0.582, 0.595, 0.608, 0.621, 0.633, 0.646, 0.659, 0.671, 0.683, 0.695, 0.707, 0.719, 0.731, 0.741, 0.750, 0.759, 0.766, 0.773, 0.778, 0.783, 0.786, 0.789, 0.791, 0.793, 0.795, 0.797, 0.798, 0.799, 0.800, 0.801, 0.801, 0.801, 0.801, 0.801, 0.800, 0.799, 0.798, 0.796, 0.794, 0.792, 0.789, 0.785, 0.782, 0.778, 0.773, 0.768, 0.763, 0.757, 0.751, 0.745, 0.738, 0.731, 0.723, 0.715, 0.707, 0.698, 0.688, 0.678, 0.668, 0.657, 0.645, 0.632, 0.620, 0.606, 0.592, 0.577, 0.562, 0.546, 0.530, 0.515, 0.499, 0.484, 0.469, 0.454, 0.439, 0.424, 0.410, 0.396, 0.381, 0.367, 0.353, 0.339, 0.324, 0.310, 0.295, 0.280, 0.265, 0.250, 0.235 };
float default_soil_vector_2[ ] = {0.253, 0.248, 0.243, 0.238, 0.232, 0.226, 0.219, 0.212, 0.205, 0.197, 0.189, 0.180, 0.171, 0.161, 0.151, 0.141, 0.130, 0.119, 0.107, 0.095, 0.082, 0.069, 0.055, 0.041, 0.026, 0.011, -.005, -.021, -.038, -.055, -.073, -.091, -.109, -.127, -.146, -.164, -.183, -.201, -.220, -.239, -.256, -.271, -.284, -.296, -.305, -.313, -.320, -.324, -.327, -.329, -.332, -.335, -.338, -.341, -.344, -.348, -.351, -.355, -.359, -.363, -.367, -.371, -.374, -.377, -.380, -.381, -.383, -.383, -.383, -.383, -.381, -.380, -.377, -.374, -.371, -.367, -.362, -.357, -.351, -.344, -.337, -.330, -.322, -.315, -.307, -.300, -.292, -.285, -.278, -.270, -.263, -.255, -.248, -.241, -.233, -.226, -.218, -.211, -.204, -.196, -.189, -.182, -.175, -.168, -.161, -.154, -.147, -.140, -.133, -.126, -.119, -.113, -.106, -.099, -.093, -.086, -.080, -.073, -.067, -.060, -.053, -.046, -.039, -.032, -.024, -.017, -.009, -.001, 0.007, 0.015, 0.023, 0.032, 0.040, 0.049, 0.058, 0.067, 0.076, 0.085, 0.095, 0.104, 0.114, 0.124, 0.134, 0.144, 0.154, 0.164, 0.175, 0.185, 0.195, 0.206, 0.216, 0.227, 0.237, 0.247, 0.258, 0.268, 0.279, 0.289, 0.300, 0.310, 0.321, 0.332, 0.342, 0.353, 0.363, 0.374, 0.385, 0.395, 0.406, 0.417, 0.428, 0.438, 0.449, 0.460, 0.471, 0.483, 0.494, 0.506, 0.518, 0.531, 0.544, 0.556, 0.570, 0.583, 0.597, 0.611, 0.625, 0.639, 0.654, 0.669, 0.684, 0.699, 0.715, 0.730, 0.745, 0.759, 0.772, 0.784, 0.796, 0.807, 0.817, 0.826, 0.834, 0.842, 0.849, 0.855, 0.860, 0.865, 0.868, 0.871, 0.873, 0.876, 0.879, 0.882, 0.885, 0.889, 0.893, 0.898, 0.903, 0.908, 0.914, 0.920, 0.926, 0.933, 0.940, 0.948, 0.956, 0.964, 0.973, 0.982, 0.990, 0.999, 1.007, 1.014, 1.022, 1.029, 1.036, 1.042, 1.048, 1.054, 1.060, 1.065, 1.070, 1.074, 1.079, 1.083, 1.086, 1.090, 1.093, 1.095, 1.098, 1.100, 1.103, 1.105, 1.108, 1.111, 1.114, 1.117, 1.120, 1.123, 1.127, 1.130, 1.134, 1.137, 1.141, 1.145, 1.149, 1.153, 1.157, 1.161, 1.165, 1.169, 1.172, 1.175, 1.177, 1.178, 1.179, 1.180, 1.180, 1.180, 1.179, 1.178, 1.176, 1.174, 1.171, 1.169, 1.166, 1.163, 1.160, 1.157, 1.153, 1.146, 1.133, 1.115, 1.089, 1.052, 1.002, 0.940, 0.874, 0.807, 0.740, 0.677, 0.624, 0.582, 0.550, 0.527, 0.510, 0.500, 0.497, 0.501, 0.507, 0.514, 0.523, 0.533, 0.544, 0.557, 0.570, 0.584, 0.598, 0.613, 0.629, 0.646, 0.663, 0.681, 0.700, 0.719, 0.738, 0.757, 0.775, 0.792, 0.809, 0.826, 0.842, 0.857, 0.872, 0.887, 0.901, 0.914, 0.927, 0.940, 0.952, 0.963, 0.974, 0.985, 0.995, 1.004, 1.013, 1.022, 1.030, 1.037, 1.044, 1.051, 1.057, 1.063, 1.068, 1.073, 1.077, 1.081, 1.084, 1.087, 1.089, 1.091, 1.092, 1.093, 1.093, 1.093, 1.092, 1.091, 1.089, 1.087, 1.085, 1.082, 1.078, 1.074, 1.070, 1.065, 1.059, 1.053, 1.047, 1.040, 1.033, 1.025, 1.017, 1.009, 0.999, 0.990, 0.980, 0.970, 0.959, 0.947, 0.935, 0.923, 0.910, 0.897, 0.883, 0.869, 0.855, 0.840, 0.824, 0.808, 0.792, 0.775, 0.757, 0.740, 0.721, 0.703, 0.683, 0.664, 0.644, 0.623, 0.602, 0.581, 0.559, 0.536, 0.513, 0.490, 0.466, 0.442, 0.417, 0.392, 0.366 };
float default_soil_vector_3[ ] = {-.455, -.412, -.369, -.327, -.286, -.245, -.204, -.164, -.124, -.085, -.046, -.008, 0.030, 0.067, 0.104, 0.140, 0.176, 0.211, 0.246, 0.280, 0.314, 0.348, 0.381, 0.413, 0.445, 0.477, 0.507, 0.538, 0.567, 0.597, 0.625, 0.653, 0.681, 0.708, 0.735, 0.761, 0.785, 0.809, 0.831, 0.852, 0.871, 0.889, 0.906, 0.921, 0.935, 0.948, 0.959, 0.969, 0.978, 0.986, 0.994, 1.000, 1.006, 1.011, 1.014, 1.017, 1.019, 1.020, 1.020, 1.019, 1.018, 1.015, 1.012, 1.008, 1.004, 0.998, 0.992, 0.985, 0.978, 0.969, 0.960, 0.950, 0.940, 0.929, 0.916, 0.904, 0.890, 0.876, 0.861, 0.845, 0.828, 0.811, 0.794, 0.775, 0.756, 0.736, 0.716, 0.695, 0.673, 0.650, 0.627, 0.605, 0.583, 0.561, 0.541, 0.520, 0.501, 0.481, 0.463, 0.444, 0.427, 0.409, 0.392, 0.375, 0.358, 0.341, 0.324, 0.307, 0.290, 0.273, 0.257, 0.240, 0.224, 0.207, 0.191, 0.174, 0.158, 0.142, 0.125, 0.109, 0.093, 0.076, 0.060, 0.043, 0.027, 0.011, -.006, -.023, -.039, -.056, -.072, -.089, -.105, -.120, -.135, -.150, -.165, -.179, -.192, -.205, -.218, -.231, -.243, -.255, -.266, -.277, -.288, -.298, -.308, -.318, -.328, -.337, -.346, -.355, -.363, -.372, -.379, -.387, -.394, -.402, -.409, -.417, -.425, -.433, -.442, -.450, -.459, -.469, -.478, -.488, -.498, -.507, -.516, -.525, -.533, -.540, -.548, -.555, -.561, -.567, -.572, -.576, -.578, -.579, -.578, -.576, -.573, -.568, -.562, -.555, -.546, -.537, -.527, -.516, -.503, -.484, -.456, -.420, -.375, -.323, -.264, -.202, -.138, -.086, -.051, -.030, -.015, -.005, -.001, 0.001, 0.002, 0.001, -.001, -.004, -.010, -.018, -.028, -.040, -.054, -.070, -.089, -.108, -.128, -.147, -.167, -.186, -.205, -.225, -.244, -.263, -.281, -.299, -.316, -.333, -.349, -.365, -.380, -.394, -.408, -.422, -.435, -.447, -.459, -.470, -.481, -.490, -.499, -.507, -.515, -.522, -.528, -.533, -.538, -.542, -.545, -.547, -.549, -.550, -.550, -.550, -.549, -.547, -.544, -.541, -.537, -.531, -.526, -.519, -.511, -.503, -.494, -.485, -.476, -.468, -.460, -.452, -.445, -.439, -.433, -.427, -.422, -.418, -.416, -.415, -.415, -.416, -.413, -.406, -.395, -.381, -.363, -.343, -.321, -.289, -.231, -.147, -.071, -.013, 0.032, 0.077, 0.122, 0.160, 0.185, 0.199, 0.203, 0.205, 0.207, 0.208, 0.209, 0.209, 0.210, 0.210, 0.209, 0.209, 0.208, 0.206, 0.204, 0.202, 0.200, 0.197, 0.194, 0.191, 0.187, 0.183, 0.178, 0.174, 0.168, 0.163, 0.157, 0.150, 0.143, 0.136, 0.129, 0.121, 0.112, 0.104, 0.095, 0.086, 0.077, 0.068, 0.060, 0.051, 0.043, 0.034, 0.026, 0.018, 0.009, 0.001, -.007, -.015, -.022, -.030, -.037, -.043, -.049, -.054, -.059, -.063, -.067, -.071, -.074, -.076, -.078, -.079, -.080, -.081, -.081, -.080, -.079, -.077, -.075, -.073, -.070, -.067, -.064, -.060, -.056, -.052, -.047, -.042, -.037, -.031, -.025, -.019, -.012, -.005, 0.003, 0.010, 0.019, 0.027, 0.036, 0.045, 0.055, 0.065, 0.076, 0.087, 0.098, 0.110, 0.122, 0.134, 0.144, 0.153, 0.160, 0.165, 0.168, 0.170, 0.171, 0.169, 0.166, 0.162, 0.157, 0.151, 0.144, 0.137, 0.128, 0.119, 0.109, 0.098, 0.086, 0.074, 0.061 }; 
float default_soil_vector_4[ ] = {0.058, 0.070, 0.081, 0.092, 0.103, 0.112, 0.122, 0.131, 0.139, 0.148, 0.155, 0.162, 0.169, 0.175, 0.181, 0.186, 0.191, 0.196, 0.200, 0.203, 0.206, 0.209, 0.211, 0.213, 0.214, 0.215, 0.215, 0.215, 0.214, 0.213, 0.212, 0.210, 0.207, 0.204, 0.201, 0.197, 0.192, 0.187, 0.182, 0.176, 0.169, 0.162, 0.154, 0.146, 0.138, 0.129, 0.119, 0.109, 0.099, 0.089, 0.079, 0.069, 0.059, 0.049, 0.039, 0.029, 0.019, 0.009, -.001, -.011, -.021, -.031, -.041, -.051, -.060, -.069, -.077, -.085, -.092, -.099, -.105, -.111, -.116, -.121, -.125, -.129, -.133, -.136, -.138, -.140, -.142, -.143, -.145, -.147, -.148, -.150, -.151, -.153, -.154, -.156, -.157, -.159, -.160, -.162, -.163, -.164, -.165, -.165, -.165, -.165, -.165, -.164, -.163, -.162, -.161, -.159, -.157, -.155, -.152, -.149, -.146, -.142, -.138, -.133, -.128, -.122, -.116, -.109, -.102, -.094, -.085, -.076, -.067, -.057, -.046, -.035, -.024, -.012, 0.000, 0.012, 0.024, 0.036, 0.049, 0.061, 0.074, 0.087, 0.100, 0.113, 0.126, 0.140, 0.153, 0.167, 0.181, 0.195, 0.210, 0.224, 0.238, 0.253, 0.267, 0.282, 0.297, 0.311, 0.326, 0.341, 0.356, 0.371, 0.385, 0.400, 0.415, 0.430, 0.446, 0.461, 0.477, 0.493, 0.510, 0.527, 0.544, 0.562, 0.579, 0.597, 0.614, 0.631, 0.648, 0.664, 0.680, 0.695, 0.710, 0.725, 0.738, 0.752, 0.764, 0.776, 0.787, 0.797, 0.805, 0.811, 0.815, 0.816, 0.816, 0.816, 0.818, 0.821, 0.820, 0.811, 0.793, 0.765, 0.720, 0.658, 0.579, 0.494, 0.408, 0.322, 0.242, 0.181, 0.141, 0.119, 0.108, 0.106, 0.113, 0.124, 0.134, 0.142, 0.149, 0.157, 0.169, 0.184, 0.204, 0.228, 0.255, 0.283, 0.313, 0.344, 0.377, 0.411, 0.446, 0.480, 0.514, 0.547, 0.579, 0.610, 0.640, 0.669, 0.698, 0.726, 0.752, 0.777, 0.801, 0.824, 0.846, 0.867, 0.886, 0.904, 0.921, 0.937, 0.952, 0.967, 0.980, 0.992, 1.004, 1.015, 1.025, 1.034, 1.042, 1.050, 1.056, 1.062, 1.067, 1.071, 1.074, 1.077, 1.078, 1.079, 1.077, 1.074, 1.069, 1.063, 1.054, 1.044, 1.033, 1.019, 1.004, 0.988, 0.971, 0.954, 0.937, 0.921, 0.906, 0.893, 0.882, 0.872, 0.865, 0.858, 0.851, 0.844, 0.838, 0.830, 0.823, 0.814, 0.805, 0.790, 0.761, 0.716, 0.655, 0.579, 0.491, 0.394, 0.287, 0.171, 0.045, -.078, -.177, -.252, -.303, -.330, -.333, -.331, -.329, -.327, -.326, -.324, -.323, -.321, -.319, -.316, -.313, -.309, -.304, -.298, -.293, -.288, -.282, -.277, -.271, -.265, -.257, -.247, -.236, -.222, -.207, -.190, -.172, -.154, -.137, -.122, -.107, -.093, -.080, -.068, -.055, -.042, -.028, -.013, 0.004, 0.021, 0.039, 0.058, 0.078, 0.095, 0.109, 0.119, 0.127, 0.130, 0.131, 0.128, 0.121, 0.112, 0.102, 0.092, 0.084, 0.075, 0.067, 0.060, 0.053, 0.047, 0.042, 0.037, 0.032, 0.028, 0.025, 0.021, 0.018, 0.015, 0.011, 0.008, 0.004, 0.001, -.003, -.006, -.010, -.014, -.017, -.021, -.025, -.028, -.032, -.036, -.040, -.044, -.047, -.051, -.054, -.056, -.059, -.060, -.062, -.063, -.063, -.063, -.063, -.063, -.062, -.060, -.058, -.056, -.053, -.050, -.047, -.045, -.042, -.040, -.039, -.037, -.035, -.034, -.033, -.033, -.032, -.032, -.032, -.032, -.033 };

#endif
