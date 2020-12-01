/* Author: Malachi Griffith
*  Date: Sept. 16 2002
*  Purpose: Given the velocity and time by a user,
*           this program calculates the distance travelled in feet.
*/

main()
{
	float time, velocity, distance;

	printf("Please enter the amount of time taken >");
	scanf("%f", &time);
	printf("Please enter your velocity >");
	scanf("%f", &velocity);	
	distance = time * velocity;
	printf("\nVelocity is %.2f ft per second\n", velocity);
	printf("time is %.1f seconds\n", time);
	printf("The distance convered is %.2f feet \n", distance);
}


