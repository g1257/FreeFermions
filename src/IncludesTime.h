#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"

typedef double RealType;
typedef std::complex<double> FieldType;
typedef std::vector<bool> LevelsType;
typedef psimag::Matrix<FieldType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,LevelsType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
