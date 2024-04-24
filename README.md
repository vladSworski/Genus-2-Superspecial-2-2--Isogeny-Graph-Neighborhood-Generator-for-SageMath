# Genus-2-Superspecial-2-2--Isogeny-Graph-Neighborhood-Generator-for-SageMath
(The following is copied directly from my PhD dissertation, "Number of 4-Cycles of the Genus 2 Superspecial Isogeny Graph".  I've edited some parts to make more sense in this context, but otherwise it's 1-1.  I'm new to GitHub and how to properly upload a project like this, so expect a LOT of changes in the near future.  The hope for at least the moment being is that this can be used to snag my files, upload them to CoCalc, and see for yourself that this works/use this as a basis for your own project.  In it's current incarnation, I have no reason to believe that this will run outside of CoCalc.  Further testing is required on my part.

This code was run in CoCalc, on SageMath version 9.8 within a Ubuntu 20.04 environment.  At the time of writing this paper, the code was in version 1.5.

\subsection*{Overview}
General calculation within the author's code is dependent on attaching/loading five `.sage' files.  In brief, they are as follows.
\begin{itemize}
\item `LGGClass.sage' This file defines classes to represent our vertices: one for Elliptic product vertices and one for Hyperelliptic Jacobian vertices.  It also maintains data in dictionaries and arrays about vertices in the graph for easy recall.
\item `LGGIsog.sage'  This file handles the calculation of isogenies between vertices.  Using the methodologies discussed in the previous sections, it can construct Richelot Isogenies.  It is also able to construct isogenies between elliptic curves, or from a product of elliptic curves to a Jacobian of Hyperelliptic curves.  It also maintains code to determine the value of $\delta$ from the previous section.
\item `LGGUtilNew.sage' This file contains a variety of functions relevant to the project as whole.
\item `VTools.sage' This file contains a variety of functions that the author considers useful in a broader context.
\item `LGGDataGen.sage' This files contains a variety of functions that allow us to calculate data on the local neighborhood graphs we produce.
\end{itemize}

The remaining files in the project are used for direct calculation in the form of `.sagews' files.
\begin{itemize}
\item `LGGScript.sagews' This file constructs and saves local neighborhood graphs to a file.
\item `LGGDataCollector.sagews' This file runs data collection tools on saved local neighborhood graphs, then saves the data for processing to file and if requested generates visualizations of the network graphs themselves.
\item `LGGMultiOutput.sagews' and `LGGDataCollectorMulti.sagews' do the same as above but process many neighborhoods in series.
\item `Type 6 Neighborhood.sagews' This file provides a walkthrough for the proof of the crab graph in the results chapter of my dissertation and the number of four cycles through the type 6 vertex.
\item `Type Sigma 1728 Neighborhood.sagews' This file provides a walkthrough for the proof of the nautilus graph in the results chapter of my dissertation and the number of four cycles through the type $\Sigma_{1728}$ vertex.
\item `Type Pi 0 1728 Neighborhood.sagews' This file provides a walkthrough for the proof of the turtle graph in the results chapter of my dissertation and the number of four cycles through the type $\Pi_{0,1728}$ vertex.
\end{itemize}

\subsection*{Constructing Neighborhoods}
The first stage of our algorithm begins by constructing local neighborhood graphs.  The `LGGScript' file or equivalent loads the necessary .sage files and packages and then the user inputs the starting curve, prime, and parameters used for limiting program run time.\par
The algorithm takes the starting curve, builds a class object for it from the `LGGClass' file and detects nearby neighbors by calculating every possible $2$-isogeny and image thereof.  It also calculates the curve's invariants for easier storage, and the corresponding automorphism group to then calculate its type.  The neighbors are added to an array for processing, and each one has the same process applied to it - adding their neighbors to the same list.  This process continues until there are no more vertices in the graph or an arbitrary limit is hit (maximum number, distance away from the starting vertex, etc.)\par
Every isogeny is calculated by calling upon the `LGGIsog' file and the appropriate method for the domain and image of the isogeny.\par
At the end of the process, the data on the neighborhoods is saved to a file using python's `pickle' package.

\subsection*{Forming Data}
We can take our saved neighborhoods, `unpickle' them, and process data from them.  This process is primarily carried out in the `LGGDataCollection' file or equivalent.\par
There are a number of things our code can do:
\begin{itemize}
\item Split data based on characteristic, so that we can analyze specific things like spinal vertices, or a specific type of vertices.
\item Calculate various types of graph theoretic matrices.
\item Calculate connectivity of subgraphs.
\item Calculate four cycles in the graph.
\item Generate a visualization of the neighborhood
\end{itemize}

Some of these tasks are worth breaking down.  Most of the functions run to accomplish these tasks are found in `LGGDataGen'.\par
One of the most important tools for splitting data is `splitSetToolSpinality' which separates spinal and non-spinal vertices into two sets.  We could also split our dataset into sets based on how far they are from the starting vertex using `splitSetToolLayer'.\par
We can calculate the adjacency matrix, diagonal matrix, and Laplace matrix for a given graph using `calculateAdjacencyMatrices'.  It is simple from here to calculate the connectivity of the data using `spineComponents' on the Laplace matrix.  There is a theorem in graph theory, that the number of components of a network graph is equal to the rank of the kernel of its Laplace matrix.  This method calculates and returns the rank.  If the number of components is anything other than 1, it isn't connected.\par
Calculating the number of four-cycles through the starting vertex is more complicated and is addressed at the end of this section.\par
We generate visualizations of neighborhoods using `plotNbhd', a method using network graph tools from python's `networkx' package.  By default, the graphs made by this method use the `circular layout'.  But `kamada kawai' can be called by including the parameter layout = 1 and `spring layout' by layout = 2.  Most of our example graphs are created using spring layout.  Networkx is not good at double edges or weights, so those are not included in our visualizations.

\subsection*{Proving Theorems}
The files used for proving the theorems in the results chapter of my dissertation are formatted as `Type \_\_\_ Neighborhood'.  These files attempt to construct the complete radius 2 neighborhood of a vertex in the graph over quadratic extensions of $\mathbb{Q}$.
These will walk you through the steps taken to construct the graph, clarifying the proofs in the results chapter of my dissertation.  Running the cells in the sage worksheets in order accomplishes this.\par
The process goes as follows.  We construct the starting vertex.  Determine if any isogenies require a field extension.  If they do, we extend our field to include that.  We calculate all isogenies and their images, picking up the outgoing weights along the way.  We calculate the invariants and automorphism group of the images and from that get their types.  We repeat this process for each of the new vertices at a distance of one away from the starting vertex.  Because we calculate invariants, we are able to determine if any vertices have common neighbors.  Because we calculate the automorphism group we are able to calculate all reverse edge weights in the graph using Smith's Edge-Weight Orbit Stabilizer Theorem.  This is enough to calculate all relevant details to the proofs.


 \subsection*{Algorithm for Detecting $4$-Cycles in the Code}
 Here is the methodology for identifying $4$-cycles in the graph in our code.\par
 `fourCyclesData' is called on the vertices and their adjacency matrix. This method  does the following:
\begin{enumerate}
\item Call the adjacency matrix `A', find its shape $(n,n).$  
\item Build a new matrix `P' with shape $(n,(n-1)(n-2)/2)$ where each column contains 0 in every row except two of them, which contain a 1.  (The 1 cannot be in the 0th row).  There is a column in the matrix for each distinct way to do this.
\item Calculate $B := A\cdot P.$  It will have shape $(n,(n-1)(n-2)/2)$.
\item For each column in $B,$ check if the the entry in the zeroth row is $2.$  If it is, check if any other position in the column has a $2.$  Record each such column.  This can be translated into the number of unweighted-undirected cycles
\item Use edge weights on this data to count the number of weighted-directed cycles.
\end{enumerate}
\begin{proposition}
There is a $4$-cycle between the nodes $0,m,z,n$ if and only if the column of B representing the edge between $m$ and $n$ contains a $2$ in the $0$'th and $z$'th rows.
\end{proposition}
\subsubsection*{An Example}
\par A four-cycle in an adjacency matrix appears as a rectangle with 1's in the corner, for example:

\[A =
\begin{bmatrix}
0&\textcolor{red}{1}&0&0&\textcolor{red}{1}\\
\textcolor{teal}{1}&0&0&\textcolor{teal}{1}&0\\
0&0&0&1&0\\
0&\textcolor{red}{1}&1&0&\textcolor{red}{1}\\
\textcolor{teal}{1}&0&0&\textcolor{teal}{1}&0
\end{bmatrix}
\]
Here we see a four cycle is present on vertices 0,1,3 and 4. The associated `P' matrix looks like this:

\[P =
\begin{bmatrix}
0&0&0&0&0&0\\
1&1&1&0&0&0\\
1&0&0&1&1&0\\
0&1&0&1&0&1\\
0&0&1&0&1&1\\
\end{bmatrix}
\]

If we take the product of these two matrices we get

\[B =
\begin{bmatrix}
1&1&2&0&1&1\\
0&1&0&1&0&1\\
0&1&0&1&0&1\\
2&1&2&1&2&1\\
0&1&0&1&0&1
\end{bmatrix}
\]

The only column that has more than a single $2$ is the result of multiplying by the input column $[0,1,0,0,1]$ and getting the image column $[2,0,0,2,0].$  We can decode this to get that a $4$-cycle through the $0$ vertex includes the $1$ vertex and $4$ vertex as direct neighbors (the position of the $1$'s in the initial vector are $1$ and $4$.) and the $3$ vertex opposite it in the cycle (the position of the $2$'s in the final vector are $0$ and $3$.)  Every such cycle we record adds to the count of unweighted-undirected $4$-cycles. From here, we can use edge weight data to uncover the exact weighted-directed number of $4$-cycles.\par
\subsubsection*{Closing Thoughts}
A four cycle going through the $0$ vertex, a.k.a the starting vertex, must have two neighbors in that cycle.  If we multiply the adjacency matrix by the vector that contains a $1$ in each of those places but $0$'s everywhere else, the $0$ position will return a $2$ as the $0$ vertex has $2$ neighbors in that vector.  But if another $2$ appears in the image, that means that there is another vertex that both of those vertices are adjacent to.  Therefore, these four vertices must form a $4$-cycle.  Note that because we want the two vertices that are adjacent to the starting vertex to be something other than the starting vertex, we construct P to have 0's in the first row. Every other combination of two vertices has its own representative column in the P matrix. The number of ways to form such representatives is $(n-1)(n-2)/2$, hence the P matrix's width.
