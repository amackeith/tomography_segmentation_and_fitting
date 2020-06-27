import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("denoise_lentils_22_select_center.py not passed enough arguments")
    exit()

else:
    fname = sys.argv[1]
    arr = np.load(fname)
    print("Loaded", fname)
    

selected_circles = []
coords = []
locked = False
finalized = False
center_x = "Not Selected Yet"
center_y = "Not Selected Yet"

slice_index = arr.shape[0] // 2
print(slice_index)

half_min_side_length = min([arr.shape[1], arr.shape[2]]) / 2.0
# the smaller of the non-verticle array sides


x = [int(half_min_side_length), ];  # x and y location on plot of approx center
y = [int(half_min_side_length), ];  # usally it is quite close to 200,200
scale = [300, ];


def get_s():
    global x, half_min_side_length
    return len(x) * half_min_side_length / 3.0 + 2


def onclick(event):
    global x, y, scale, finalized, locked
    if finalized:
        return
    
    s = get_s()  # scale for this particular instance
    
    if event.button == 1 and event.inaxes != None:
        # guess a location (will place a cirlce )
        finalized = False
        print("Guess: x: %s y: %s" % (event.xdata, event.ydata))
        print('you have guessed this many times:', len(x) + 1)
        # left click adds a circle
        print(event)
        print(event.button)
        scale.append(np.pi * (s ** 2))
        x.append(event.xdata)
        y.append(event.ydata)
        # clear frame
        plt.clf()
        plt.imshow(arr[slice_index], cmap='gray', origin='lower')
        plt.xlabel("e1")
        plt.ylabel("e2")
        plt.title("left click to add circles, (recomend using at least 5)\n" +
                  "right click to remove last circle or approve mean")
        plt.scatter(x, y, scale, facecolors='none', edgecolors='r')
        # inform matplotlib of the new data
        plt.draw()  # redraw
    
    if event.button == 3:
        # right clicking gives you the option to say your are done,
        # you want to remove a circle, or that it was a mistake
        done = {'d', 'done'}
        no = {'no', 'n'}
        remove = {'r', 'remove'}
        sys.stdout.write("\nPress done to lock in selected centers\n" +
                         "remove to remove previouse circle\n" +
                         "no to add more circles\n")
        choice = input().lower()
        if choice in remove:
            print("Remove Last Click")
            if len(x) >= 1:
                x = x[:-1]
                y = y[:-1]
                scale = scale[:-1]
            
            # clear frame
            plt.clf()
            plt.imshow(arr[slice_index], cmap='gray', origin='lower')
            plt.xlabel("e1")
            plt.ylabel("e2")
            plt.title(
                "left click to add circles, (recomend using at least 5)\n" +
                " right click to remove last circle or approve mean")
            plt.scatter(x, y, scale, facecolors='none', edgecolors='r')
            # inform matplotlib of the new data
            plt.draw()  # redraw
        elif choice in no:
            plt.clf()
            plt.imshow(arr[slice_index], cmap='gray', origin='lower')
            plt.xlabel("e1")
            plt.ylabel("e2")
            plt.title("left click to add circles, (recomend using at least 5)" +
                      "\n right click to remove last circle or approve mean")
            plt.scatter(x, y, scale, facecolors='none', edgecolors='r')
            # inform matplotlib of the new data
            plt.draw()  # redraw
        elif choice in done:
            finalized = True
            
            best_x = [np.mean(x), np.mean(x)]
            best_y = [np.mean(y), np.mean(y)]
            plt.clf()
            plt.imshow(arr[slice_index], cmap='gray', origin='lower')
            plt.xlabel("e1")
            plt.ylabel("e2")
            plt.title("To approve, close the plot, else right click\
             and select no")
            plt.scatter(best_x, best_y, [20, 30], facecolors='none',
                        edgecolors='r');  # inform matplotlib of the new data
            plt.draw()  # redraw
        
        else:
            sys.stdout.write("Please respond with 'done', 'remove', or 'no'\n")
        
        # right click removes the previouse circle


def motion_notify(event):
    # this places the purple ring and lets you place it where it makes sense
    global x, y, scale, locked, finalized
    if locked or finalized:
        return
    s = get_s()  # scale for this particular instance
    s = np.pi * (s ** 2)
    # clear frame and plot new stuff this time
    # including a hypothetical circle that may or may not be added
    plt.clf()
    plt.imshow(arr[slice_index], cmap='gray', origin='lower')
    plt.xlabel("e1")
    plt.ylabel("e2")
    plt.title(
        "left click to add circles, \
        (recomend using at least 5)\n \
         right click to remove last circle or approve mean")
    plt.scatter([event.xdata], [event.ydata], s, facecolors='none'
                , edgecolors='purple')  # adds candidate cirlce
    
    if len(scale)!=0:
        plt.scatter(x, y, scale, facecolors='none',
                    edgecolors='r')  # inform matplotlib of the new data
    plt.draw()  # redraw


def handle_close(event):
    global finalized, locked
    if finalized == False:
        sys.stdout.write("ERROR: Program not\
         properly exited bad mean returned\n")
    
    elif len(x) == 0:
        sys.stdout.write("ERROR: No Circles have been selected,\
         please select the center of the rings\n")
    
    else:
        sys.stdout.write("Good Exit, mean over %s points e1: %s ,e2: %s\n" %
                         (len(x), np.mean(x), np.mean(y)))
        locked = True
        out = np.array([np.mean(x), np.mean(y)])
        np.save("center_of_ring.npy", out)


fig, ax = plt.subplots()
fig.canvas.mpl_connect('button_release_event', onclick)
fig.canvas.mpl_connect('motion_notify_event', motion_notify)
fig.canvas.mpl_connect('close_event', handle_close)

plt.title("Line up circles with ring artifact to\n\
            determine it's center for de noising left Click to begin")
plt.xlabel("e1")
plt.ylabel("e2")
plt.show()
# plt.draw()
