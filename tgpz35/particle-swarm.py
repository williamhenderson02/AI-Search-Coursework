############
############ ALTHOUGH I GIVE YOU THIS TEMPLATE PROGRAM WITH THE NAME 'skeleton.py', 
############ YOU CAN RENAME IT TO ANYTHING YOU LIKE. HOWEVER, FOR THE PURPOSES OF 
############ THE EXPLANATION IN THESE COMMENTS, I ASSUME THAT THIS PROGRAM IS STILL 
############ CALLED 'skeleton.py'.
############
############ IF YOU WISH TO IMPORT STANDARD MODULES, YOU CAN ADD THEM AFTER THOSE BELOW.
############ NOTE THAT YOU ARE NOT ALLOWED TO IMPORT ANY NON-STANDARD MODULES! TO SEE
############ THE STANDARD MODULES, TAKE A LOOK IN 'validate_before_handin.py'.
############
############ DO NOT INCLUDE ANY COMMENTS ON A LINE WHERE YOU IMPORT A MODULE.
############

import os
import sys
import time
import random
import math
import decimal

############ START OF SECTOR 1 (IGNORE THIS COMMENT)
############
############ NOW PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS.
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############ BY 'DO NOT TOUCH' I REALLY MEAN THIS. EVEN CHANGING THE SYNTAX, BY
############ ADDING SPACES OR COMMENTS OR LINE RETURNS AND SO ON, COULD MEAN THAT
############ CODES MIGHT NOT RUN WHEN I RUN THEM!
############

def read_file_into_string(input_file, ord_range):
    the_file = open(input_file, 'r') 
    current_char = the_file.read(1) 
    file_string = ""
    length = len(ord_range)
    while current_char != "":
        i = 0
        while i < length:
            if ord(current_char) >= ord_range[i][0] and ord(current_char) <= ord_range[i][1]:
                file_string = file_string + current_char
                i = length
            else:
                i = i + 1
        current_char = the_file.read(1)
    the_file.close()
    return file_string

def remove_all_spaces(the_string):
    length = len(the_string)
    new_string = ""
    for i in range(length):
        if the_string[i] != " ":
            new_string = new_string + the_string[i]
    return new_string

def integerize(the_string):
    length = len(the_string)
    stripped_string = "0"
    for i in range(0, length):
        if ord(the_string[i]) >= 48 and ord(the_string[i]) <= 57:
            stripped_string = stripped_string + the_string[i]
    resulting_int = int(stripped_string)
    return resulting_int

def convert_to_list_of_int(the_string):
    list_of_integers = []
    location = 0
    finished = False
    while finished == False:
        found_comma = the_string.find(',', location)
        if found_comma == -1:
            finished = True
        else:
            list_of_integers.append(integerize(the_string[location:found_comma]))
            location = found_comma + 1
            if the_string[location:location + 5] == "NOTE=":
                finished = True
    return list_of_integers

def build_distance_matrix(num_cities, distances, city_format):
    dist_matrix = []
    i = 0
    if city_format == "full":
        for j in range(num_cities):
            row = []
            for k in range(0, num_cities):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    elif city_format == "upper_tri":
        for j in range(0, num_cities):
            row = []
            for k in range(j):
                row.append(0)
            for k in range(num_cities - j):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    else:
        for j in range(0, num_cities):
            row = []
            for k in range(j + 1):
                row.append(0)
            for k in range(0, num_cities - (j + 1)):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    if city_format == "upper_tri" or city_format == "strict_upper_tri":
        for i in range(0, num_cities):
            for j in range(0, num_cities):
                if i > j:
                    dist_matrix[i][j] = dist_matrix[j][i]
    return dist_matrix

def read_in_algorithm_codes_and_tariffs(alg_codes_file):
    flag = "good"
    code_dictionary = {}   
    tariff_dictionary = {}  
    if not os.path.exists(alg_codes_file):
        flag = "not_exist"  
        return code_dictionary, tariff_dictionary, flag
    ord_range = [[32, 126]]
    file_string = read_file_into_string(alg_codes_file, ord_range)  
    location = 0
    EOF = False
    list_of_items = []  
    while EOF == False: 
        found_comma = file_string.find(",", location)
        if found_comma == -1:
            EOF = True
            sandwich = file_string[location:]
        else:
            sandwich = file_string[location:found_comma]
            location = found_comma + 1
        list_of_items.append(sandwich)
    third_length = int(len(list_of_items)/3)
    for i in range(third_length):
        code_dictionary[list_of_items[3 * i]] = list_of_items[3 * i + 1]
        tariff_dictionary[list_of_items[3 * i]] = int(list_of_items[3 * i + 2])
    return code_dictionary, tariff_dictionary, flag

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY!
############
############ THE RESERVED VARIABLE 'input_file' IS THE CITY FILE UNDER CONSIDERATION.
############
############ IT CAN BE SUPPLIED BY SETTING THE VARIABLE BELOW OR VIA A COMMAND-LINE
############ EXECUTION OF THE FORM 'python skeleton.py city_file.txt'. WHEN SUPPLYING
############ THE CITY FILE VIA A COMMAND-LINE EXECUTION, ANY ASSIGNMENT OF THE VARIABLE
############ 'input_file' IN THE LINE BELOW iS SUPPRESSED.
############
############ IT IS ASSUMED THAT THIS PROGRAM 'skeleton.py' SITS IN A FOLDER THE NAME OF
############ WHICH IS YOUR USER-NAME, E.G., 'abcd12', WHICH IN TURN SITS IN ANOTHER
############ FOLDER. IN THIS OTHER FOLDER IS THE FOLDER 'city-files' AND NO MATTER HOW
############ THE NAME OF THE CITY FILE IS SUPPLIED TO THIS PROGRAM, IT IS ASSUMED THAT 
############ THE CITY FILE IS IN THE FOLDER 'city-files'.
############
############ END OF SECTOR 1 (IGNORE THIS COMMENT)

input_file = "AISearchfile012.txt"

############ START OF SECTOR 2 (IGNORE THIS COMMENT)
############
############ PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS STARTING
############ 'HAVE YOU TOUCHED ...'
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############

if len(sys.argv) > 1:
    input_file = sys.argv[1]

############ END OF SECTOR 2 (IGNORE THIS COMMENT)
path_for_city_files = "../city-files"
############ START OF SECTOR 3 (IGNORE THIS COMMENT)
    
if os.path.isfile(path_for_city_files + "/" + input_file):
    ord_range = [[32, 126]]
    file_string = read_file_into_string(path_for_city_files + "/" + input_file, ord_range)
    file_string = remove_all_spaces(file_string)
    print("I have found and read the input file " + input_file + ":")
else:
    print("*** error: The city file " + input_file + " does not exist in the city-file folder.")
    sys.exit()

location = file_string.find("SIZE=")
if location == -1:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()
    
comma = file_string.find(",", location)
if comma == -1:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()
    
num_cities_as_string = file_string[location + 5:comma]
num_cities = integerize(num_cities_as_string)
print("   the number of cities is stored in 'num_cities' and is " + str(num_cities))

comma = comma + 1
stripped_file_string = file_string[comma:]
distances = convert_to_list_of_int(stripped_file_string)

counted_distances = len(distances)
if counted_distances == num_cities * num_cities:
    city_format = "full"
elif counted_distances == (num_cities * (num_cities + 1))/2:
    city_format = "upper_tri"
elif counted_distances == (num_cities * (num_cities - 1))/2:
    city_format = "strict_upper_tri"
else:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()

dist_matrix = build_distance_matrix(num_cities, distances, city_format)
print("   the distance matrix 'dist_matrix' has been built.")

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY!
############
############ YOU NOW HAVE THE NUMBER OF CITIES STORED IN THE INTEGER VARIABLE 'num_cities'
############ AND THE TWO_DIMENSIONAL MATRIX 'dist_matrix' HOLDS THE INTEGER CITY-TO-CITY 
############ DISTANCES SO THAT 'dist_matrix[i][j]' IS THE DISTANCE FROM CITY 'i' TO CITY 'j'.
############ BOTH 'num_cities' AND 'dist_matrix' ARE RESERVED VARIABLES AND SHOULD FEED
############ INTO YOUR IMPLEMENTATIONS.
############
############ THERE NOW FOLLOWS CODE THAT READS THE ALGORITHM CODES AND TARIFFS FROM
############ THE TEXT-FILE 'alg_codes_and_tariffs.txt' INTO THE RESERVED DICTIONARIES
############ 'code_dictionary' AND 'tariff_dictionary'. DO NOT AMEND THIS CODE!
############ THE TEXT FILE 'alg_codes_and_tariffs.txt' SHOULD BE IN THE SAME FOLDER AS
############ THE FOLDER 'city-files' AND THE FOLDER WHOSE NAME IS YOUR USER-NAME.
############
############ PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS STARTING
############ 'HAVE YOU TOUCHED ...'
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############
############ END OF SECTOR 3 (IGNORE THIS COMMENT)

############ START OF SECTOR 4 (IGNORE THIS COMMENT)
path_for_alg_codes_and_tariffs = "../alg_codes_and_tariffs.txt"
############ END OF SECTOR 4 (IGNORE THIS COMMENT)

############ START OF SECTOR 5 (IGNORE THIS COMMENT)
code_dictionary, tariff_dictionary, flag = read_in_algorithm_codes_and_tariffs(path_for_alg_codes_and_tariffs)

if flag != "good":
    print("*** error: The text file 'alg_codes_and_tariffs.txt' does not exist.")
    sys.exit()

print("The codes and tariffs have been read from 'alg_codes_and_tariffs.txt':")

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY! SORRY TO GO ON ABOUT THIS BUT YOU NEED TO BE 
############ AWARE OF THIS FACT!
############
############ YOU NOW NEED TO SUPPLY SOME PARAMETERS.
############
############ THE RESERVED STRING VARIABLE 'my_user_name' SHOULD BE SET AT YOUR
############ USER-NAME, E.G., "abcd12"
############
############ END OF SECTOR 5 (IGNORE THIS COMMENT)

my_user_name = "tgp35"

############ START OF SECTOR 6 (IGNORE THIS COMMENT)
############
############ YOU CAN SUPPLY, IF YOU WANT, YOUR FULL NAME. THIS IS NOT USED AT ALL BUT SERVES AS
############ AN EXTRA CHECK THAT THIS FILE BELONGS TO YOU. IF YOU DO NOT WANT TO SUPPLY YOUR
############ NAME THEN EITHER SET THE STRING VARIABLES 'my_first_name' AND 'my_last_name' AT 
############ SOMETHING LIKE "Mickey" AND "Mouse" OR AS THE EMPTY STRING (AS THEY ARE NOW;
############ BUT PLEASE ENSURE THAT THE RESERVED VARIABLES 'my_first_name' AND 'my_last_name'
############ ARE SET AT SOMETHING).
############
############ END OF SECTOR 6 (IGNORE THIS COMMENT)

my_first_name = "William"
my_last_name = "Henderson"

############ START OF SECTOR 7 (IGNORE THIS COMMENT)
############
############ YOU NEED TO SUPPLY THE ALGORITHM CODE IN THE RESERVED STRING VARIABLE 'algorithm_code'
############ FOR THE ALGORITHM YOU ARE IMPLEMENTING. IT NEEDS TO BE A LEGAL CODE FROM THE TEXT-FILE
############ 'alg_codes_and_tariffs.txt' (READ THIS FILE TO SEE THE CODES).
############
############ END OF SECTOR 7 (IGNORE THIS COMMENT)

algorithm_code = "PS"

############ START OF SECTOR 8 (IGNORE THIS COMMENT)
############
############ PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS STARTING
############ 'HAVE YOU TOUCHED ...'
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############

if not algorithm_code in code_dictionary:
    print("*** error: the algorithm code " + algorithm_code + " is illegal")
    sys.exit()
print("   your algorithm code is legal and is " + algorithm_code + " -" + code_dictionary[algorithm_code] + ".")

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY! SORRY TO GO ON ABOUT THIS BUT YOU NEED TO BE 
############ AWARE OF THIS FACT!
############
############ YOU CAN ADD A NOTE THAT WILL BE ADDED AT THE END OF THE RESULTING TOUR FILE IF YOU LIKE,
############ E.G., "in my basic greedy search, I broke ties by always visiting the first 
############ city found" BY USING THE RESERVED STRING VARIABLE 'added_note' OR LEAVE IT EMPTY
############ IF YOU WISH. THIS HAS NO EFFECT ON MARKS BUT HELPS YOU TO REMEMBER THINGS ABOUT
############ YOUR TOUR THAT YOU MIGHT BE INTERESTED IN LATER.
############
############ END OF SECTOR 8 (IGNORE THIS COMMENT)

added_note = ""

############
############ NOW YOUR CODE SHOULD BEGIN.
############

def pso(max_it, N, delta):

    def initialise_positions(N):

        particles = []

        start_city = random.randint(0, num_cities - 1)

        for i in range(N):

            tour = [start_city]
            unvisited_cities = [j for j in range(num_cities) if j != start_city]

            for j in range(num_cities-1):

                city = random.choice(unvisited_cities)
                tour.append(city)
                unvisited_cities.remove(city)

            particles.append(tour)

        return particles

    def initialise_velocities(N):

        velocities = []

        for i in range(N):

            velocity = []

            #potentially add a max number of times a city can be swapped
            for j in range(num_cities - 1):
                #index 1 so first city is never swapped
                swap_index = random.randint(1, num_cities - 1)  # choose one of the cities in the swap

                #last city has 2 ways to get swapped and not just 1
                #downside is penultimate city now has 3 ways to be swapped so more likeley to be swapped than other cities
                #better than all cities having twice the chance as last city
                if swap_index == num_cities - 1:
                
                    swap = swap_index - 1, swap_index

                else: 

                    swap = swap_index, swap_index + 1  # represents the swap of two consecutive cities

                velocity.append(swap)


            velocities.append(velocity)

        return velocities

    def get_min_tour(tours):

        best_length = math.inf
        best_tour = []

        for particle in tours:

            tour_length = 0

            for i in range(0, len(particle)):

                tour_length += dist_matrix[particle[i-1]][particle[i]]

            if tour_length < best_length:
                best_length = tour_length
                best_tour = particle
        
        return best_tour

    def get_metric_distance(particle_a, particle_b):
        
        swaps = 0
        swap_tuple = ()
        linear_order = particle_b.copy()
        sorting = particle_a.copy()
        velocity  = []

        is_sorted = False

        while not is_sorted:

            is_sorted = True

            for i in range(1, len(particle_a) - 1):

                num_1 = sorting[i]
                num_2 = sorting[i+1]

                index_1 = linear_order.index(num_1)
                index_2 = linear_order.index(num_2)

                if index_1 > index_2:

                    swaps += 1

                    swap_tuple = i, i + 1
                    velocity.append(swap_tuple)

                    sorting[i] = num_2
                    sorting[i+1] = num_1

                    is_sorted = False

        return swaps, velocity

    def get_neighbourhood(particles,particle):

        particle_index = particles.index(particle)
        potential_neighbours = [particles[i] for i in range(len(particles)) if i != particle_index]
        neighbourhood = []

        for neighbour in potential_neighbours:

            distance, velocity = get_metric_distance(particle, neighbour)

            if distance <= delta:

                neighbourhood.append(neighbour)
            
        return neighbourhood

    def get_n_best(neighbourhood):

        lengths = []

        for particle in neighbourhood:

            tour_length = 0

            for i in range(0, len(particle)):

                tour_length += dist_matrix[particle[i-1]][particle[i]]

            lengths.append(tour_length)

        n_best = min(lengths)

        return get_n_best

    def transform_particle_position(particle, velocitity):

        particle_transformed = particle.copy()

        for i in velocity:

            index = i[0]
            particle_transformed[index], particle_transformed[index + 1] = particle_transformed[index + 1], particle_transformed[index]

        return particle_transformed

    def compose_particle_velocity(particle, velocity, p_hat, n_best):

        particle_swaps, particle_contribution = get_metric_distance(particle, p_hat)

        if n_best != math.inf:

            neighbourhood_swaps, neighbourhood_contribution = get_metric_distance(particle, n_best)

        else:

            neighbourhood_swaps = 0
            neighbourhood_contribution = []

        return particle_swaps, particle_contribution, neighbourhood_swaps, neighbourhood_contribution

    def inertia_function(t, w_start, w_end):

        w = w_start - ((w_start - w_end) * t) / max_it

        return w

    def multiply_velocity(weight, velocity):

        length = len(velocity)

        if weight == 1:

            return velocity

        elif weight < 1:

            nums_to_take = int((length) * (weight))
            velocity = velocity[0:nums_to_take]

            return velocity

        else:

            decimal_weight = decimal.Decimal(str(weight))
            integer_part = int(decimal_weight)
            decimal_part = float(decimal_weight - integer_part)
            decimal_velocity = velocity.copy()
            integer_velocity = velocity.copy()

            for i in range(1,integer_part):

                if len(integer_velocity) != 0:

                    for j in range(len(integer_velocity)):

                        velocity.append(integer_velocity[j])

            if decimal_part != 0:

                nums_to_take = int((length) * (decimal_part))
                decimal_velocity = decimal_velocity[0:nums_to_take]

            else:

                decimal_velocity = None

            if decimal_velocity != None:

                for i in range(len(decimal_velocity)):

                    velocity.append(decimal_velocity[i])

            return velocity

    def calc_next_velocity(inertia, cognitive_factor, social_factor):

        next_velocity = []

        if len(inertia) != 0:

            for i in range(len(inertia)):

                next_velocity.append(inertia[i])

        if len(cognitive_factor) != 0:

            for i in range(len(cognitive_factor)):

                next_velocity.append(cognitive_factor[i])

        if len(social_factor) != 0:

            for i in range(len(social_factor)):

                next_velocity.append(social_factor[i])

        return next_velocity

    particles = initialise_positions(N)
    p_hats = particles.copy()
    velocities = initialise_velocities(N)
    p_best = get_min_tour(p_hats)

    tour_length = 0
    
    for i in range(0, len(p_best)):

        tour_length += dist_matrix[p_best[i-1]][p_best[i]]

    print("start length", tour_length)

    t = 0
    w_start = 0.9
    w_end = 0.4
    alpha = 0.75
    beta = 2.75

    while t < max_it:
        #print("t", t)

        possible_bests = [p_best]

        for particle in particles:

            index = particles.index(particle)
            velocity = velocities[index]
            neighbourhood = get_neighbourhood(particles,particle)
            p_hat = p_hats[index]

            if len(neighbourhood) != 0:

                n_best = get_min_tour(neighbourhood)

            else:

                n_best = math.inf

            next_position = transform_particle_position(particle,velocity)

            particle_swaps, particle_contribution, neighbourhood_swaps, neighbourhood_contribution = compose_particle_velocity(particle, velocity, p_hat, n_best)

            inertia_weight = inertia_function(t, w_start, w_end)

            if len(particle_contribution) != 0:
                
                epsilon = random.randint(0,len(particle_contribution) - 1)

            else:

                epsilon = 0

            if len(neighbourhood_contribution) != 0:

                epsilon_prime = random.randint(0, len(neighbourhood_contribution) - 1)

            else:

                epsilon_prime = 0

            if particle_swaps != 0:

                particle_contribution.pop(epsilon)
            
            if neighbourhood_swaps != 0:

                neighbourhood_contribution.pop(epsilon_prime)

            inertia = multiply_velocity(inertia_weight, velocity)

            cognitive_factor = multiply_velocity(alpha, particle_contribution)

            social_factor = multiply_velocity(beta, neighbourhood_contribution)
        
            next_velocity = calc_next_velocity(inertia, cognitive_factor, social_factor)

            compare_tours = []
            compare_tours.append(next_position)
            compare_tours.append(p_hat)

            next_p_hat = get_min_tour(compare_tours)
            possible_bests.append(next_p_hat)

            particles[index] = next_position
            velocities[index] = next_velocity
            p_hats[index] = next_p_hat

        p_best = get_min_tour(possible_bests)

        t += 1

    tour_length = 0
    
    for i in range(0, len(p_best)):

        tour_length += dist_matrix[p_best[i-1]][p_best[i]]

    print("end length", tour_length)

    return p_best

pso(1000,10,20)



############ START OF SECTOR 9 (IGNORE THIS COMMENT)
############
############ YOUR CODE SHOULD NOW BE COMPLETE AND WHEN EXECUTION OF THIS PROGRAM 'skeleton.py'
############ REACHES THIS POINT, YOU SHOULD HAVE COMPUTED A TOUR IN THE RESERVED LIST VARIABLE 'tour', 
############ WHICH HOLDS A LIST OF THE INTEGERS FROM {0, 1, ..., 'num_cities' - 1} SO THAT EVERY INTEGER
############ APPEARS EXACTLY ONCE, AND YOU SHOULD ALSO HOLD THE LENGTH OF THIS TOUR IN THE RESERVED
############ INTEGER VARIABLE 'tour_length'.
############
############ YOUR TOUR WILL BE PACKAGED IN A TOUR FILE OF THE APPROPRIATE FORMAT AND THIS TOUR FILE'S,
############ NAME WILL BE A MIX OF THE NAME OF THE CITY FILE, THE NAME OF THIS PROGRAM AND THE
############ CURRENT DATE AND TIME. SO, EVERY SUCCESSFUL EXECUTION GIVES A TOUR FILE WITH A UNIQUE
############ NAME AND YOU CAN RENAME THE ONES YOU WANT TO KEEP LATER.
############
############ DO NOT TOUCH OR ALTER THE CODE BELOW THIS POINT! YOU HAVE BEEN WARNED!
############

"""flag = "good"
length = len(tour)
for i in range(0, length):
    if isinstance(tour[i], int) == False:
        flag = "bad"
    else:
        tour[i] = int(tour[i])
if flag == "bad":
    print("*** error: Your tour contains non-integer values.")
    sys.exit()
if isinstance(tour_length, int) == False:
    print("*** error: The tour-length is a non-integer value.")
    sys.exit()
tour_length = int(tour_length)
if len(tour) != num_cities:
    print("*** error: The tour does not consist of " + str(num_cities) + " cities as there are, in fact, " + str(len(tour)) + ".")
    sys.exit()
flag = "good"
for i in range(0, num_cities):
    if not i in tour:
        flag = "bad"
if flag == "bad":
    print("*** error: Your tour has illegal or repeated city names.")
    sys.exit()
check_tour_length = 0
for i in range(0, num_cities - 1):
    check_tour_length = check_tour_length + dist_matrix[tour[i]][tour[i + 1]]
check_tour_length = check_tour_length + dist_matrix[tour[num_cities - 1]][tour[0]]
if tour_length != check_tour_length:
    flag = print("*** error: The length of your tour is not " + str(tour_length) + "; it is actually " + str(check_tour_length) + ".")
    sys.exit()
print("You, user " + my_user_name + ", have successfully built a tour of length " + str(tour_length) + "!")

local_time = time.asctime(time.localtime(time.time()))
output_file_time = local_time[4:7] + local_time[8:10] + local_time[11:13] + local_time[14:16] + local_time[17:19]
output_file_time = output_file_time.replace(" ", "0")
script_name = os.path.basename(sys.argv[0])
if len(sys.argv) > 2:
    output_file_time = sys.argv[2]
output_file_name = script_name[0:len(script_name) - 3] + "_" + input_file[0:len(input_file) - 4] + "_" + output_file_time + ".txt"

f = open(output_file_name,'w')
f.write("USER = " + my_user_name + " (" + my_first_name + " " + my_last_name + "),\n")
f.write("ALGORITHM CODE = " + algorithm_code + ", NAME OF CITY-FILE = " + input_file + ",\n")
f.write("SIZE = " + str(num_cities) + ", TOUR LENGTH = " + str(tour_length) + ",\n")
f.write(str(tour[0]))
for i in range(1,num_cities):
    f.write("," + str(tour[i]))
f.write(",\nNOTE = " + added_note)
f.close()
print("I have successfully written your tour to the tour file:\n   " + output_file_name + ".")"""

############ END OF SECTOR 9 (IGNORE THIS COMMENT)
    
    











    


