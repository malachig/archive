/* Author: Malachi Griffith
*  Date: Sept. 16 2002
*  Purpose: Given the velocity and time this program calculates distance
*/

main()
{
	float time, velocity, distance;

	time = 0.3;
	velocity = 24.0;
	distance = time * velocity;
	printf("Velocity is %.2f ft per second\n", velocity);
	printf("time is %.1f seconds\n", time);
	printf("The distance convered is %.2f feet \n", distance);
}


