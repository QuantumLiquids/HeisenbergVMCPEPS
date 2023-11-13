% system size: 8x8; method: fixed step length stochastic gradient
e0 = [-39.34742960, -38.33357654, -39.22839946, -39.32205511, -39.34805267, -39.35908875, -39.36911554, -39.37332720, -39.40016801, -39.40365964, -39.41497824, -39.42242562, -39.43092539, -39.43544149, -39.43956988, -39.44272606, -39.44876759, -39.45627259, -39.45968109, -39.45767312, -39.46351429, -39.47033675, -39.46896048, -39.47333709, -39.47755187, -39.48085076, -39.48419375, -39.48352382, -39.48350219, -39.49165746];
e0_persite = e0/64;

e0_ex = -0.619040;

semilogy((e0_persite -e0_ex)/abs(e0_ex), 'o','LineWidth', 2); hold on;