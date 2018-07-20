#include "vtkPointMatcher.h"
#include "vtkPointDistanceMatrix.h"
#include "vtkCombinatoricGenerator.h"
#include <vtkDoubleArray.h>
#include <vtkLandmarkTransform.h>
#include <vtkSmartPointer.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMath.h>

#define RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR VTK_DOUBLE_MAX
#define MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH 3
#define MAXIMUM_NUMBER_OF_POINTS_NEEDED_FOR_DETERMINISTIC_MATCH 6

//----------------------------------------------------------------------------
vtkStandardNewMacro( vtkPointMatcher );

//------------------------------------------------------------------------------
vtkPointMatcher::vtkPointMatcher()
{
  this->InputPointList1 = NULL;
  this->InputPointList2 = NULL;
  this->MaximumDifferenceInNumberOfPoints = 2;
  this->TolerableRootMeanSquareDistanceErrorMm = 10.0;
  this->ComputedRootMeanSquareDistanceErrorMm = RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR;
  this->AmbiguityThresholdDistanceMm = 5.0;
  this->MatchingAmbiguous = false;
  // outputs are never null
  this->OutputPointList1 = vtkSmartPointer< vtkPoints >::New();
  this->OutputPointList2 = vtkSmartPointer< vtkPoints >::New();

  // timestamps for input and output are the same, initially
  this->Modified();
  this->OutputChangedTime.Modified();
}

//------------------------------------------------------------------------------
vtkPointMatcher::~vtkPointMatcher()
{
}

//------------------------------------------------------------------------------
void vtkPointMatcher::PrintSelf( std::ostream &os, vtkIndent indent )
{
  Superclass::PrintSelf( os, indent );
  
  os << indent << "MaximumDifferenceInNumberOfPoints: " << this->MaximumDifferenceInNumberOfPoints << std::endl;
  os << indent << "TolerableRootMeanSquareDistanceErrorMm: " << this->TolerableRootMeanSquareDistanceErrorMm << std::endl;
  os << indent << "ComputedRootMeanSquareDistanceErrorMm: " << this->ComputedRootMeanSquareDistanceErrorMm << std::endl;
  os << indent << "IsMatchingWithinTolerance: " << this->IsMatchingWithinTolerance() << std::endl;
  os << indent << "IsMatchingAmbiguous" << this->IsMatchingAmbiguous() << std::endl;
  os << indent << "UpdateNeeded: " << this->UpdateNeeded() << std::endl;
}

//------------------------------------------------------------------------------
// INPUT MUTATORS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void vtkPointMatcher::SetInputPointList1( vtkPoints* points )
{
  vtkSetObjectBodyMacro( InputPointList1, vtkPoints, points );
  // mean distance error has not been computed yet. Set to maximum possible value for now
  this->ComputedRootMeanSquareDistanceErrorMm = RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR;
}

//------------------------------------------------------------------------------
void vtkPointMatcher::SetInputPointList2( vtkPoints* points )
{
  vtkSetObjectBodyMacro( InputPointList2, vtkPoints, points );
  // mean distance error has not been computed yet. Set to maximum possible value for now
  this->ComputedRootMeanSquareDistanceErrorMm = RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR;
}

//------------------------------------------------------------------------------
void vtkPointMatcher::SetMaximumDifferenceInNumberOfPoints( unsigned int numberOfPoints )
{
  this->MaximumDifferenceInNumberOfPoints = numberOfPoints;
  this->Modified();
  // mean distance error has not been computed yet. Set to maximum possible value for now
  this->ComputedRootMeanSquareDistanceErrorMm = RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR;
}

//------------------------------------------------------------------------------
void vtkPointMatcher::SetTolerableRootMeanSquareDistanceErrorMm( double errorMm )
{
  this->TolerableRootMeanSquareDistanceErrorMm = errorMm;
  this->Modified();
  // mean distance error has not been computed yet. Set to maximum possible value for now
  this->ComputedRootMeanSquareDistanceErrorMm = RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR;
}

//------------------------------------------------------------------------------
void vtkPointMatcher::SetAmbiguityThresholdDistanceMm( double thresholdMm )
{
  this->AmbiguityThresholdDistanceMm = thresholdMm;
  this->Modified();
  // mean distance error has not been computed yet. Set to maximum possible value for now
  this->ComputedRootMeanSquareDistanceErrorMm = RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR;
}

//------------------------------------------------------------------------------
// OUTPUT ACCESSORS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
vtkPoints* vtkPointMatcher::GetOutputPointList1()
{
  if ( this->UpdateNeeded() )
  {
    this->Update();
  }

  return this->OutputPointList1;
}

//------------------------------------------------------------------------------
vtkPoints* vtkPointMatcher::GetOutputPointList2()
{
  if ( this->UpdateNeeded() )
  {
    this->Update();
  }

  return this->OutputPointList2;
}

//------------------------------------------------------------------------------
double vtkPointMatcher::GetComputedRootMeanSquareDistanceErrorMm()
{
  if ( this->UpdateNeeded() )
  {
    this->Update();
  }

  return this->ComputedRootMeanSquareDistanceErrorMm;
}

//------------------------------------------------------------------------------
bool vtkPointMatcher::IsMatchingWithinTolerance()
{
  if ( this->UpdateNeeded() )
  {
    this->Update();
  }
  
  return ( this->ComputedRootMeanSquareDistanceErrorMm <= this->TolerableRootMeanSquareDistanceErrorMm );
}

//------------------------------------------------------------------------------
bool vtkPointMatcher::IsMatchingAmbiguous()
{
  if ( this->UpdateNeeded() )
  {
    this->Update();
  }

  return this->MatchingAmbiguous;
}

//------------------------------------------------------------------------------
// LOGIC
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void vtkPointMatcher::Update()
{
  if ( !this->UpdateNeeded() )
  {
    return;
  }
  
  this->ComputedRootMeanSquareDistanceErrorMm = RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR;
  this->MatchingAmbiguous = false;
  this->OutputPointList1->Reset();
  this->OutputPointList2->Reset();

  // check number of points, make sure point lists are already similar/same sizes
  if ( this->InputPointList1 == NULL )
  {
    vtkWarningMacro( "Input point list 1 is null. Cannot update." );
    return;
  }
  int pointList1Size = this->InputPointList1->GetNumberOfPoints();

  if ( this->InputPointList2 == NULL )
  {
    vtkWarningMacro( "Input point list 2 is null. Cannot update." );
    return;
  }
  int pointList2Size = this->InputPointList2->GetNumberOfPoints();
  
  unsigned int differenceInPointListSizes = abs( pointList1Size - pointList2Size );
  if ( differenceInPointListSizes > this->MaximumDifferenceInNumberOfPoints )
  {
    // no matching to do... just output the first N points of both lists
    int smallestNumberOfPoints = vtkMath::Min( pointList1Size, pointList2Size );
    vtkPointMatcher::CopyFirstNPoints( this->InputPointList1, this->OutputPointList1, smallestNumberOfPoints );
    vtkPointMatcher::CopyFirstNPoints( this->InputPointList2, this->OutputPointList2, smallestNumberOfPoints );
    this->OutputChangedTime.Modified();
    return;
  }

  int numberOfPointsInList1 = this->InputPointList1->GetNumberOfPoints();
  int numberOfPointsInList2 = this->InputPointList2->GetNumberOfPoints();
  int smallerPointListSize = vtkMath::Min( numberOfPointsInList1, numberOfPointsInList2 );
  if ( numberOfPointsInList1 < MAXIMUM_NUMBER_OF_POINTS_NEEDED_FOR_DETERMINISTIC_MATCH &&
       numberOfPointsInList2 < MAXIMUM_NUMBER_OF_POINTS_NEEDED_FOR_DETERMINISTIC_MATCH )
  {
    // point set is small enough for brute force approach
    int minimumSubsetSize = vtkMath::Max( ( smallerPointListSize - this->MaximumDifferenceInNumberOfPoints ), ( unsigned int )MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH );
    int maximumSubsetSize = smallerPointListSize;
    vtkPointMatcher::UpdateBestMatchingForSubsetsOfPoints( minimumSubsetSize, maximumSubsetSize,
                                                            this->InputPointList1, this->InputPointList2,
                                                            this->AmbiguityThresholdDistanceMm, this->MatchingAmbiguous,
                                                            this->ComputedRootMeanSquareDistanceErrorMm, this->TolerableRootMeanSquareDistanceErrorMm,
                                                            this->OutputPointList1, this->OutputPointList2 );
  }
  else
  {
    // Point set is large so brute force is unsuitable

    // Sort points according to the 'uniqueness' of their geometry relative to other points
    vtkSmartPointer< vtkPoints > inputPointList1SortedByUniqueness = vtkSmartPointer< vtkPoints >::New();
    vtkPointMatcher::ReorderPointsAccordingToUniqueGeometry( this->InputPointList1, inputPointList1SortedByUniqueness );
    vtkSmartPointer< vtkPoints > inputPointList2SortedByUniqueness = vtkSmartPointer< vtkPoints >::New();
    vtkPointMatcher::ReorderPointsAccordingToUniqueGeometry( this->InputPointList2, inputPointList2SortedByUniqueness );

    // Take the first N points that are most 'unique'
    int numberOfPointsToUseForInitialRegistration = vtkMath::Min( smallerPointListSize, MAXIMUM_NUMBER_OF_POINTS_NEEDED_FOR_DETERMINISTIC_MATCH );
    vtkSmartPointer< vtkPoints > inputPointList1Reduced = vtkSmartPointer< vtkPoints >::New();
    vtkPointMatcher::CopyFirstNPoints( inputPointList1SortedByUniqueness, inputPointList1Reduced, numberOfPointsToUseForInitialRegistration );
    vtkSmartPointer< vtkPoints > inputPointList2Reduced = vtkSmartPointer< vtkPoints >::New();
    vtkPointMatcher::CopyFirstNPoints( inputPointList2SortedByUniqueness, inputPointList2Reduced, numberOfPointsToUseForInitialRegistration );

    // Compute correspondence between those points
    int numberOfPointsInReducedList1 = inputPointList1Reduced->GetNumberOfPoints();
    int numberOfPointsInReducedList2 = inputPointList2Reduced->GetNumberOfPoints();
    int smallerPointListSize = vtkMath::Min( numberOfPointsInReducedList1, numberOfPointsInReducedList2 );
    int minimumSubsetSize = vtkMath::Max( ( smallerPointListSize - this->MaximumDifferenceInNumberOfPoints ), ( unsigned int )MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH );
    int maximumSubsetSize = smallerPointListSize;
    vtkSmartPointer< vtkPoints > matchedPointList1 = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkPoints > matchedPointList2 = vtkSmartPointer< vtkPoints >::New();
    vtkPointMatcher::UpdateBestMatchingForSubsetsOfPoints( minimumSubsetSize, maximumSubsetSize,
                                                           inputPointList1Reduced, inputPointList2Reduced,
                                                           this->AmbiguityThresholdDistanceMm, this->MatchingAmbiguous,
                                                           this->ComputedRootMeanSquareDistanceErrorMm, this->TolerableRootMeanSquareDistanceErrorMm,
                                                           matchedPointList1, matchedPointList2 );
    
    // Compute initial registration based on this correspondence
    vtkSmartPointer< vtkLandmarkTransform > registrationList1ToList2Transform = vtkSmartPointer< vtkLandmarkTransform >::New();
    registrationList1ToList2Transform->SetSourceLandmarks( matchedPointList1 );
    registrationList1ToList2Transform->SetTargetLandmarks( matchedPointList2 );
    registrationList1ToList2Transform->SetModeToRigidBody();
    matchedPointList2->Modified();
    registrationList1ToList2Transform->Update();

    vtkSmartPointer< vtkPolyData > inputPointList1PolyData = vtkSmartPointer< vtkPolyData >::New();
    inputPointList1PolyData->SetPoints( this->InputPointList1 );
    
    vtkSmartPointer< vtkTransformPolyDataFilter > registrationList1ToList2TransformFilter = vtkSmartPointer< vtkTransformPolyDataFilter >::New();
    registrationList1ToList2TransformFilter->SetTransform( registrationList1ToList2Transform );
    registrationList1ToList2TransformFilter->SetInputData( inputPointList1PolyData );
    registrationList1ToList2TransformFilter->Update();

    // iterative closest point
    int numberOfIterationsRemaining = 100;
    vtkPoints* registeredPointList1; // reused
    vtkSmartPointer< vtkPoints > nearestPointsList2 = vtkSmartPointer< vtkPoints >::New(); // reused
    vtkSmartPointer< vtkPointDistanceMatrix > pointDistanceMatrix = vtkSmartPointer< vtkPointDistanceMatrix >::New(); // reused
    vtkSmartPointer< vtkDoubleArray > pointToPointDistances = vtkSmartPointer< vtkDoubleArray >::New(); // reused
    while ( numberOfIterationsRemaining > 0 )
    {
      // perform the registration
      registrationList1ToList2TransformFilter->Update();

      // get the transformed list 1
      vtkPolyData* registeredPointList1PolyData = vtkPolyData::SafeDownCast( registrationList1ToList2TransformFilter->GetOutput() );
      if ( registeredPointList1PolyData == NULL )
      {
        vtkGenericWarningMacro( "Matched point list poly data is null." );
        return;
      }

      registeredPointList1 = vtkPoints::SafeDownCast( registeredPointList1PolyData->GetPoints() );
      if ( registeredPointList1 == NULL )
      {
        vtkGenericWarningMacro( "Matched point list is null." );
        return;
      }

      // find nearest points from list 2
      nearestPointsList2->Reset();
      for ( int pointList1Index = 0; pointList1Index < numberOfPointsInList1; pointList1Index++ )
      {
        int nearestPointList2Index = 0;
        double distanceSquaredToNearestPoint = VTK_DOUBLE_MAX;
        double transformedPointFromList1[ 3 ];
        registeredPointList1->GetPoint( pointList1Index, transformedPointFromList1 );
        for ( int pointList2Index = 0; pointList2Index < numberOfPointsInList2; pointList2Index++ )
        {
          double pointFromList2[ 3 ];
          this->InputPointList2->GetPoint( pointList2Index, pointFromList2 );
          double distanceSquaredToPointFromList2 = vtkMath::Distance2BetweenPoints( transformedPointFromList1, pointFromList2 );
          if ( distanceSquaredToPointFromList2 < distanceSquaredToNearestPoint )
          {
            nearestPointList2Index = pointList2Index;
            distanceSquaredToNearestPoint = distanceSquaredToPointFromList2;
          }
        }
        double nearestPoint[ 3 ];
        this->InputPointList2->GetPoint( nearestPointList2Index, nearestPoint );
        nearestPointsList2->InsertNextPoint( nearestPoint );
      }

      // check how good this registration is
      pointDistanceMatrix->SetPointList1( registeredPointList1 );
      pointDistanceMatrix->SetPointList2( nearestPointsList2 );
      pointDistanceMatrix->Update();

      pointToPointDistances->Reset();
      pointDistanceMatrix->GetDistances( pointToPointDistances );
      int numberOfDistances = pointToPointDistances->GetNumberOfTuples();
      if ( numberOfDistances == 0 )
      {
        vtkGenericWarningMacro( "There are no distances to determine quality of registration." );
        return;
      }
      
      double sumDistancesSquared = 0.0;
      for ( int distanceIndex = 0; distanceIndex < numberOfDistances; distanceIndex++ )
      {
        double distance = pointToPointDistances->GetComponent( distanceIndex, 0 );
        double distanceSquared = distance * distance;
        sumDistancesSquared += distanceSquared;
      }
      this->ComputedRootMeanSquareDistanceErrorMm = sqrt( sumDistancesSquared / numberOfDistances );

      if ( this->ComputedRootMeanSquareDistanceErrorMm < this->TolerableRootMeanSquareDistanceErrorMm )
      {
        break;
      }
      
      // update the registration parameters
      registrationList1ToList2Transform->SetSourceLandmarks( registeredPointList1 );
      registrationList1ToList2Transform->SetTargetLandmarks( nearestPointsList2 );
      nearestPointsList2->Modified();

      numberOfIterationsRemaining -= 1;
    }
    this->OutputPointList1->DeepCopy( registeredPointList1 );
    this->OutputPointList2->DeepCopy( nearestPointsList2 );
  }

  this->OutputChangedTime.Modified();
}

//------------------------------------------------------------------------------
void vtkPointMatcher::UpdateBestMatchingForSubsetsOfPoints( int minimumSubsetSize,
                                                            int maximumSubsetSize,
                                                            vtkPoints* unmatchedPointList1,
                                                            vtkPoints* unmatchedPointList2,
                                                            double ambiguityThresholdDistanceMm,
                                                            bool& matchingAmbiguous, 
                                                            double& computedRootMeanSquareDistanceErrorMm,
                                                            double tolerableRootMeanSquareDistanceErrorMm,
                                                            vtkPoints* outputMatchedPointList1,
                                                            vtkPoints* outputMatchedPointList2 )
{
  // lots of error checking
  if ( maximumSubsetSize < MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH )
  {
    vtkGenericWarningMacro( "Maximum subset size " << minimumSubsetSize << " must be at least " << MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH << ". Setting to " << MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH << "." );
    maximumSubsetSize = MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH;
  }

  if ( minimumSubsetSize < MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH )
  {
    vtkGenericWarningMacro( "Minimum subset size " << minimumSubsetSize << " must be at least " << MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH << ". Setting to " << MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH << "." );
    minimumSubsetSize = MINIMUM_NUMBER_OF_POINTS_NEEDED_TO_MATCH;
  }

  if ( maximumSubsetSize < minimumSubsetSize )
  {
    vtkGenericWarningMacro( "Maximum subset size " << maximumSubsetSize << " must be greater than or equal to minimum subset size " << minimumSubsetSize << ". Setting to " << minimumSubsetSize << "." );
    maximumSubsetSize = minimumSubsetSize;
  }

  if ( unmatchedPointList1 == NULL )
  {
    vtkGenericWarningMacro( "Unmatched point list 1 is null." );
    return;
  }

  int numberOfPointsInUnmatched1 = unmatchedPointList1->GetNumberOfPoints();
  if ( numberOfPointsInUnmatched1 < maximumSubsetSize )
  {
    vtkGenericWarningMacro( "Maximum subset size is " << maximumSubsetSize << " but the unmatched point list 1 has only " << numberOfPointsInUnmatched1 << " points." );
    return;
  }

  if ( unmatchedPointList2 == NULL )
  {
    vtkGenericWarningMacro( "Unmatched point list 2 is null." );
    return;
  }

  int numberOfPointsInUnmatched2 = unmatchedPointList2->GetNumberOfPoints();
  if ( numberOfPointsInUnmatched2 < maximumSubsetSize )
  {
    vtkGenericWarningMacro( "Maximum subset size is " << maximumSubsetSize << " but the unmatched point list 2 has only " << numberOfPointsInUnmatched2 << " points." );
    return;
  }

  if ( outputMatchedPointList1 == NULL )
  {
    vtkGenericWarningMacro( "Output matched point list 1 is null." );
    return;
  }

  if ( outputMatchedPointList2 == NULL )
  {
    vtkGenericWarningMacro( "Output matched point list 2 is null." );
    return;
  }

  for ( int subsetSize = maximumSubsetSize; subsetSize >= minimumSubsetSize; subsetSize-- )
  {
    vtkPointMatcher::UpdateBestMatchingForNSizedSubsetsOfPoints( subsetSize,
                                                                 unmatchedPointList1, unmatchedPointList2,
                                                                 ambiguityThresholdDistanceMm, matchingAmbiguous,
                                                                 computedRootMeanSquareDistanceErrorMm,
                                                                 outputMatchedPointList1, outputMatchedPointList2 );
    if ( computedRootMeanSquareDistanceErrorMm <= tolerableRootMeanSquareDistanceErrorMm )
    {
      // suitable solution has been found, no need to continue searching
      break;
    }
  }
}

//------------------------------------------------------------------------------
void vtkPointMatcher::UpdateBestMatchingForNSizedSubsetsOfPoints(
  int subsetSize,
  vtkPoints* pointList1,
  vtkPoints* pointList2,
  double ambiguityThresholdDistanceMm,
  bool& matchingAmbiguous,
  double& computedRootMeanSquareDistanceErrorMm,
  vtkPoints* outputMatchedPointList1,
  vtkPoints* outputMatchedPointList2 )
{
  if ( pointList1 == NULL )
  {
    vtkGenericWarningMacro( "Point list 1 is null." );
    return;
  }

  int pointList1NumberOfPoints = pointList1->GetNumberOfPoints();
  if ( pointList1NumberOfPoints < subsetSize )
  {
    vtkGenericWarningMacro( "Looking for subsets of " << subsetSize << " points, but list 1 has only " << pointList1NumberOfPoints << " points." );
    return;
  }

  if ( pointList2 == NULL )
  {
    vtkGenericWarningMacro( "Point list 2 is null." );
    return;
  }

  int pointList2NumberOfPoints = pointList2->GetNumberOfPoints();
  if ( pointList2NumberOfPoints < subsetSize )
  {
    vtkGenericWarningMacro( "Looking for subsets of " << subsetSize << " points, but list 2 has only " << pointList2NumberOfPoints << " points." );
    return;
  }

  // generate sets of indices for all possible combinations of both input sets
  vtkSmartPointer< vtkCombinatoricGenerator > pointList1CombinationGenerator = vtkSmartPointer< vtkCombinatoricGenerator >::New();
  pointList1CombinationGenerator->SetCombinatoricToCombination();
  pointList1CombinationGenerator->SetSubsetSize( subsetSize );
  pointList1CombinationGenerator->SetNumberOfInputSets( 1 );
  for ( int pointIndex = 0; pointIndex < pointList1NumberOfPoints; pointIndex++ )
  {
    pointList1CombinationGenerator->AddInputElement( 0, pointIndex );
  }
  pointList1CombinationGenerator->Update();
  std::vector< std::vector< int > > pointList1CombinationIndices = pointList1CombinationGenerator->GetOutputSets();

  vtkSmartPointer< vtkCombinatoricGenerator > pointList2CombinationGenerator = vtkSmartPointer< vtkCombinatoricGenerator >::New();
  pointList2CombinationGenerator->SetCombinatoricToCombination();
  pointList2CombinationGenerator->SetSubsetSize( subsetSize );
  pointList2CombinationGenerator->SetNumberOfInputSets( 1 );
  for ( int pointIndex = 0; pointIndex < pointList2NumberOfPoints; pointIndex++ )
  {
    pointList2CombinationGenerator->AddInputElement( 0, pointIndex );
  }
  pointList2CombinationGenerator->Update();
  std::vector< std::vector< int > > pointList2CombinationIndices = pointList2CombinationGenerator->GetOutputSets();

  // these will store the actual combinations of points themselves (not indices)
  vtkSmartPointer< vtkPoints > pointList1Combination = vtkSmartPointer< vtkPoints >::New();
  pointList1Combination->SetNumberOfPoints( subsetSize );
  vtkSmartPointer< vtkPoints > pointList2Combination = vtkSmartPointer< vtkPoints >::New();
  pointList2Combination->SetNumberOfPoints( subsetSize );

  // iterate over all combinations of both input point sets
  for ( unsigned int pointList1CombinationIndex = 0; pointList1CombinationIndex < pointList1CombinationIndices.size(); pointList1CombinationIndex++ )
  {
    for ( unsigned int pointList2CombinationIndex = 0; pointList2CombinationIndex < pointList2CombinationIndices.size(); pointList2CombinationIndex++ )
    {
      // store appropriate contents in the pointList1Combination and pointList2Combination variables
      for ( vtkIdType pointIndex = 0; pointIndex < subsetSize; pointIndex++ )
      {
        vtkIdType point1Index = ( vtkIdType ) pointList1CombinationIndices[ pointList1CombinationIndex ][ pointIndex ];
        double* pointFromList1 = pointList1->GetPoint( point1Index );
        pointList1Combination->SetPoint( pointIndex, pointFromList1 );

        vtkIdType point2Index = ( vtkIdType ) pointList2CombinationIndices[ pointList2CombinationIndex ][ pointIndex ];
        double* pointFromList2 = pointList2->GetPoint( point2Index );
        pointList2Combination->SetPoint( pointIndex, pointFromList2 );
      }
      // finally see how good this particular combination is
      vtkPointMatcher::UpdateBestMatchingForSubsetOfPoints( pointList1Combination, pointList2Combination,
                                                            ambiguityThresholdDistanceMm, matchingAmbiguous,
                                                            computedRootMeanSquareDistanceErrorMm,
                                                            outputMatchedPointList1, outputMatchedPointList2 );
    }
  }
}

//------------------------------------------------------------------------------
// point pair matching will be based on the distances between each pair of ordered points.
// we have two input point lists. We want to reorder the second list such that the 
// point-to-point distances are as close as possible to those in the first.
// We will permute over all possibilities (and only ever keep the best result.)
void vtkPointMatcher::UpdateBestMatchingForSubsetOfPoints(
  vtkPoints* pointSubset1, 
  vtkPoints* pointSubset2,
  double ambiguityThresholdDistanceMm,
  bool& matchingAmbiguous,
  double& computedRootMeanSquareDistanceErrorMm,
  vtkPoints* outputMatchedPointList1,
  vtkPoints* outputMatchedPointList2 )
{
  // error checking
  if ( pointSubset1 == NULL )
  {
    vtkGenericWarningMacro( "Reference point set is NULL. This is a coding error. Please report." );
    return;
  }

  if ( pointSubset2 == NULL )
  {
    vtkGenericWarningMacro( "Compare point set is NULL. This is a coding error. Please report." );
    return;
  }

  int numberOfPoints = pointSubset1->GetNumberOfPoints(); // sizes of both lists should be identical
  if ( numberOfPoints != pointSubset2->GetNumberOfPoints() )
  {
    vtkGenericWarningMacro( "Point sets are of different sizes. This is a coding error. Please report." );
    return;
  }

  // Point distance matrix
  vtkSmartPointer< vtkPointDistanceMatrix > pointSubset1DistanceMatrix = vtkSmartPointer< vtkPointDistanceMatrix >::New();
  pointSubset1DistanceMatrix->SetPointList1( pointSubset1 );
  pointSubset1DistanceMatrix->SetPointList2( pointSubset1 ); // distances to itself
  pointSubset1DistanceMatrix->Update();

  // compute the permutations, store them in a vtkIntArray.
  vtkSmartPointer< vtkCombinatoricGenerator > combinatoricGenerator = vtkSmartPointer< vtkCombinatoricGenerator >::New();
  combinatoricGenerator->SetCombinatoricToPermutation();
  combinatoricGenerator->SetSubsetSize( numberOfPoints );
  combinatoricGenerator->SetNumberOfInputSets( 1 );
  for ( int pointIndex = 0; pointIndex < numberOfPoints; pointIndex++ )
  {
    combinatoricGenerator->AddInputElement( 0, pointIndex );
  }
  combinatoricGenerator->Update();
  std::vector< std::vector< int > > pointSubset2IndexPermutations = combinatoricGenerator->GetOutputSets();

  // iterate over all permutations - look for the most 'suitable'
  // point matching that gives distances most similar to the reference

  // create + allocate the permuted compare list once outside
  // the loop to avoid allocation/deallocation time costs.
  vtkSmartPointer< vtkPoints > permutedPointSubset2 = vtkSmartPointer< vtkPoints >::New();
  permutedPointSubset2->DeepCopy( pointSubset2 ); // fill it with placeholder data, same size as compareList
  vtkSmartPointer< vtkPointDistanceMatrix > permutedPointSubset2DistanceMatrix = vtkSmartPointer< vtkPointDistanceMatrix >::New();
  int numberOfPermutations = pointSubset2IndexPermutations.size();
  for ( int permutationIndex = 0; permutationIndex < numberOfPermutations; permutationIndex++ )
  {
    // fill permutedPointSubset2 with points from pointSubset2,
    // in the order indicate by the permuted indices
    for ( int pointIndex = 0; pointIndex < numberOfPoints; pointIndex++ )
    {
      int permutedPointIndex = pointSubset2IndexPermutations[ permutationIndex][ pointIndex ];
      double* permutedPoint = pointSubset2->GetPoint( permutedPointIndex );
      permutedPointSubset2->SetPoint( pointIndex, permutedPoint );
    }
    permutedPointSubset2DistanceMatrix->SetPointList1( permutedPointSubset2 );
    permutedPointSubset2DistanceMatrix->SetPointList2( permutedPointSubset2 );
    permutedPointSubset2DistanceMatrix->Update();
    double rootMeanSquareDistanceErrorMm = vtkPointMatcher::ComputeRootMeanSquareistanceErrors( pointSubset1DistanceMatrix, permutedPointSubset2DistanceMatrix );

    // case analysis for setting MatchingAmbiguous:
    // let ComputedRootMeanSquareDistanceErrorMm store the distance error for the *best* matching
    // let rootMeanSquareDistanceErrorMm store the distance error for the *current* matching
    // use AmbiguityThresholdDistanceMm and MatchingAmbiguous as described in the header file
    // 1
    // rootMeanSquareDistanceErrorMm is better (lower) than ComputedRootMeanSquareDistanceErrorMm,
    // but by less than AmbiguityThresholdDistanceMm
    // Result => MatchingAmbiguous should be set to true
    //   - trivial justification
    // 2
    // rootMeanSquareDistanceErrorMm is worse (higher) than ComputedRootMeanSquareDistanceErrorMm,
    // but by less than AmbiguityThresholdDistanceMm
    // Result => MatchingAmbiguous should be set to true
    //   - trivial justification
    // 3
    // rootMeanSquareDistanceErrorMm is better (lower) than ComputedRootMeanSquareDistanceErrorMm,
    // but by more than AmbiguityThresholdDistanceMm
    // Result => MatchingAmbiguous should be set to false.
    //   - If there was a _previous_ rootMeanSquareDistanceErrorMm within AmbiguityThresholdDistanceMm,
    //     that would have become ComputedRootMeanSquareDistanceErrorMm. Therefore there have not been
    //     any _previous_ rootMeanSquareDistanceErrorMm within AmbiguityThresholdDistanceMm.
    //   - If _later_ there is a rootMeanSquareDistanceErrorMm within AmbiguityThresholdDistanceMm,
    //     then MatchingAmbiguous will be set to true by either case 1 or case 2.
    //     Cases 1 and 2 will always catch an ambiguous matching, because the search is exhaustive.
    // 4
    // rootMeanSquareDistanceErrorMm is worse (higher) than ComputedRootMeanSquareDistanceErrorMm,
    // but by more than AmbiguityThresholdDistanceMm
    // Result => do nothing
    //   - This result does not matter. It *cannot* be within AmbiguityThresholdDistanceMm
    //     of the best (final) ComputedRootMeanSquareDistanceErrorMm

    // flag ambiguous if suitability within some threshold of the best result so far
    double differenceComparedToBestMm = computedRootMeanSquareDistanceErrorMm - rootMeanSquareDistanceErrorMm;
    bool withinAmbiguityThresholdDistanceMm = ( fabs( differenceComparedToBestMm ) <= ambiguityThresholdDistanceMm );

    // is this the best matching?
    if ( rootMeanSquareDistanceErrorMm < computedRootMeanSquareDistanceErrorMm )
    {
      if ( withinAmbiguityThresholdDistanceMm )
      {
        // case 1, described above
        matchingAmbiguous = true;
      }
      else
      {
        // case 3, described above
        matchingAmbiguous = false;
      }
      computedRootMeanSquareDistanceErrorMm = rootMeanSquareDistanceErrorMm;
      outputMatchedPointList1->DeepCopy( pointSubset1 );
      outputMatchedPointList2->DeepCopy( permutedPointSubset2 );
    }
    else if ( withinAmbiguityThresholdDistanceMm )
    {
      // case 2, described above
      matchingAmbiguous = true;
    }
  }
}

//------------------------------------------------------------------------------
double vtkPointMatcher::ComputeRootMeanSquareistanceErrors( vtkPointDistanceMatrix* distanceMatrix1, vtkPointDistanceMatrix* distanceMatrix2 )
{
  if ( distanceMatrix1 == NULL || distanceMatrix2 == NULL )
  {
    vtkGenericWarningMacro( "One of the input distance matrices is null. Cannot compute similarity. Returning default value " << RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR << "." );
    return RESET_VALUE_COMPUTED_ROOT_MEAN_DISTANCE_ERROR;
  }

  vtkSmartPointer< vtkDoubleArray > distanceErrorMatrix = vtkSmartPointer< vtkDoubleArray >::New();
  vtkPointDistanceMatrix::ComputePairWiseDifferences( distanceMatrix2, distanceMatrix1, distanceErrorMatrix );

  double sumOfSquaredDistanceErrors = 0;
  int numberOfColumns = distanceErrorMatrix->GetNumberOfTuples();
  int numberOfRows = distanceErrorMatrix->GetNumberOfComponents();
  for ( int columnIndex = 0; columnIndex < numberOfColumns; columnIndex++ )
  {
    for ( int rowIndex = 0; rowIndex < numberOfRows; rowIndex++ )
    {
      double currentDistanceError = distanceErrorMatrix->GetComponent( columnIndex, rowIndex );
      sumOfSquaredDistanceErrors += ( currentDistanceError * currentDistanceError );
    }
  }

  int numberOfDistances = numberOfColumns * numberOfRows;
  double meanOfSquaredDistanceErrors = sumOfSquaredDistanceErrors / numberOfDistances;
  double rootMeanSquareistanceErrors = sqrt( meanOfSquaredDistanceErrors );

  return rootMeanSquareistanceErrors;
}

//------------------------------------------------------------------------------
void vtkPointMatcher::ReorderPointsAccordingToUniqueGeometry( vtkPoints* inputUnsortedPointList, vtkPoints* outputSortedPointList )
{
  if ( inputUnsortedPointList == NULL )
  {
    vtkGenericWarningMacro( "Input point list is null." );
    return;
  }

  if ( outputSortedPointList == NULL )
  {
    vtkGenericWarningMacro( "Output sorted point list is null." );
    return;
  }

  // Figure out which points are the most 'unique'
  vtkSmartPointer< vtkDoubleArray > pointUniquenesses = vtkSmartPointer< vtkDoubleArray >::New();
  vtkPointMatcher::ComputeUniquenessForPointList( inputUnsortedPointList, pointUniquenesses );

  // sanity check
  int numberOfPoints = inputUnsortedPointList->GetNumberOfPoints();
  int numberOfUniquenesses = pointUniquenesses->GetNumberOfTuples();
  if ( numberOfUniquenesses != numberOfPoints )
  {
    vtkGenericWarningMacro( "Number of points " << numberOfPoints << " does not match number of uniquenesses " << numberOfUniquenesses << ". This is a programming error, please report it." );
    return;
  }

  // Order the list in increasing order (smallest point uniqueness is the most 'unique')
  outputSortedPointList->DeepCopy( inputUnsortedPointList );
  for ( int currentPointIndex = 0; currentPointIndex < numberOfPoints; currentPointIndex++ )
  {
    double currentUniqueness = pointUniquenesses->GetComponent( currentPointIndex, 0 );
    for ( int otherPointIndex = currentPointIndex + 1; otherPointIndex < numberOfPoints; otherPointIndex++ )
    {
      double otherUniqueness = pointUniquenesses->GetComponent( otherPointIndex, 0 );
      if ( otherUniqueness < currentUniqueness )
      {
        // swap uniquenesses
        pointUniquenesses->SetComponent( otherPointIndex, 0, currentUniqueness );
        pointUniquenesses->SetComponent( currentPointIndex, 0, otherUniqueness );
        // swap points
        double currentPoint[ 3 ];
        outputSortedPointList->GetPoint( currentPointIndex, currentPoint );
        double otherPoint[ 3 ];
        outputSortedPointList->GetPoint( otherPointIndex, otherPoint );
        outputSortedPointList->SetPoint( currentPointIndex, otherPoint );
        outputSortedPointList->SetPoint( otherPointIndex, currentPoint );
      }
    }
  }
}

//------------------------------------------------------------------------------
void vtkPointMatcher::ComputeUniquenessForPointList( vtkPoints* pointList, vtkDoubleArray* pointListUniquenesses )
{
  if ( pointList == NULL )
  {
    vtkGenericWarningMacro( "Point list is null." );
    return;
  }

  if ( pointListUniquenesses == NULL )
  {
    vtkGenericWarningMacro( "Point uniqueness array is null." );
    return;
  }

  vtkSmartPointer< vtkPointDistanceMatrix > pointDistanceMatrix = vtkSmartPointer< vtkPointDistanceMatrix >::New();
  pointDistanceMatrix->SetPointList1( pointList );
  pointDistanceMatrix->SetPointList2( pointList ); // distances to self
  pointDistanceMatrix->Update();

  vtkSmartPointer< vtkDoubleArray > allDistancesArray = vtkSmartPointer< vtkDoubleArray >::New();
  pointDistanceMatrix->GetDistances( allDistancesArray );
  double maximumDistance = pointDistanceMatrix->GetMaximumDistance();

  pointListUniquenesses->Reset();
  int numberOfPoints = pointList->GetNumberOfPoints();
  for ( int pointIndex = 0; pointIndex < numberOfPoints; pointIndex++ )
  {
    // heuristic measure for point uniqueness:
    // sum of point-to-point distance uniquenesses
    double sumOfDistanceUniquenesses = 0;
    for ( int otherPointIndex = 0; otherPointIndex < numberOfPoints; otherPointIndex++ )
    {
      double currentDistance = pointDistanceMatrix->GetDistance( pointIndex, otherPointIndex );
      sumOfDistanceUniquenesses += vtkPointMatcher::ComputeUniquenessForDistance( currentDistance, maximumDistance, allDistancesArray );
    }
    double pointUniqueness = sumOfDistanceUniquenesses;
    pointListUniquenesses->InsertNextTuple1( pointUniqueness );
  }
}


//------------------------------------------------------------------------------
double vtkPointMatcher::ComputeUniquenessForDistance( double distance, double maximumDistance, vtkDoubleArray* allDistancesArray )
{
  if ( allDistancesArray == NULL )
  {
    vtkGenericWarningMacro( "Distances array is null" );
    return 0.0;
  }

  if ( maximumDistance == 0 )
  {
    vtkGenericWarningMacro( "Maximum distance is zero. Setting to 1 to avoid division by zero." );
    maximumDistance = 1.0;
  }

  double distanceUniqueness = 0.0;
  int numberOfDistances = allDistancesArray->GetNumberOfTuples();
  for ( int distanceIndex = 0; distanceIndex < numberOfDistances; distanceIndex++ )
  {
    double otherDistance = allDistancesArray->GetComponent( distanceIndex, 0 );
    double heuristicMeasure = abs( distance - otherDistance ) / maximumDistance; // will be bounded between 0..1
    distanceUniqueness += heuristicMeasure;
  }
  return distanceUniqueness;
}

//------------------------------------------------------------------------------
void vtkPointMatcher::CopyFirstNPoints( vtkPoints* inputList, vtkPoints* outputList, int n )
{
  if ( inputList == NULL )
  {
    vtkGenericWarningMacro( "Input list is null." );
    return;
  }

  int numberOfInputPoints = inputList->GetNumberOfPoints();
  if ( n > numberOfInputPoints )
  {
    vtkGenericWarningMacro( "Call to copy " << n << " points, but there are only " << numberOfInputPoints << " points in the list." );
    return;
  }

  if ( outputList == NULL )
  {
    vtkGenericWarningMacro( "Input list is null." );
    return;
  }

  outputList->Reset();

  for ( int pointIndex = 0; pointIndex < n; pointIndex++ )
  {
    double point1[ 3 ];
    inputList->GetPoint( pointIndex, point1 );
    outputList->InsertNextPoint( point1 );
  }
}

//------------------------------------------------------------------------------
bool vtkPointMatcher::UpdateNeeded()
{
  return ( this->GetMTime() > this->OutputChangedTime );
}

